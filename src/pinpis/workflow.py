from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from .config import require
from .commands import run_command, ensure_dir, require_executable
import subprocess

def setup_directories(config: dict) -> None:
    ensure_dir(require(config, "project.output_dir"))
    ensure_dir(require(config, "project.work_dir"))
    ensure_dir(require(config, "project.tmp_dir"))


def check_external_tools(config: dict) -> None:
    external = config.get("external", {})

    for key in ["bcftools", "samtools", "tabix", "bgzip", "yn00"]:
        exe = external.get(key, key)
        require_executable(exe)

    # GNU parallel is useful but not strictly required because some scripts have fallback behavior.
    parallel = external.get("parallel", "parallel")
    if not Path("/usr/bin/parallel").exists():
        # Do not fail here; just warn.
        print(f"[pinpis] Warning: GNU parallel may not be available as '{parallel}'.")


def check_input_files(config: dict) -> None:
    input_paths = [
        require(config, "inputs.vcf"),
        require(config, "inputs.gff"),
        require(config, "inputs.reference_genome"),
        require(config, "inputs.samples_id_file"),
        require(config, "inputs.sample_years_file"),
        require(config, "inputs.gene_list"),
    ]

    for path in input_paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"Input file not found: {path}")

    for group in require(config, "sampling.groups"):
        path = group["population_file"]
        if not Path(path).exists():
            raise FileNotFoundError(f"Population file not found for group {group['name']}: {path}")


def run_sampling(config: dict) -> None:
    scripts_dir = Path(require(config, "legacy.scripts_dir"))
    python = config.get("external", {}).get("python", "python3")
    multysampler = scripts_dir / "multysampler.py"

    number_of_samples = str(require(config, "sampling.number_of_samples"))
    samples_id_file = require(config, "inputs.samples_id_file")
    sample_years_file = require(config, "inputs.sample_years_file")

    for group in require(config, "sampling.groups"):
        run_command([
            python,
            str(multysampler),
            group["population_file"],
            number_of_samples,
            samples_id_file,
            sample_years_file,
            group["output_dir"],
        ])

def read_gene_ids(config: dict) -> list:
    gene_list = Path(require(config, "inputs.gene_list"))
    genes = []

    with gene_list.open() as handle:
        for line in handle:
            gene = line.strip()
            if gene:
                genes.append(gene)

    return genes


def get_gene_info_from_gff(gene: str, gff_file: str, use_gene_name: bool) -> tuple:
    """
    Return:
        gene_id, genomic_region

    genomic_region format:
        chr:start-end
    """
    pattern = f"Name={gene}" if use_gene_name else f"ID={gene}"

    with open(gff_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]

            if feature_type == "gene" and pattern in attributes:
                chrom = fields[0]
                start = fields[3]
                end = fields[4]

                if use_gene_name:
                    gene_id = None
                    for item in attributes.split(";"):
                        if item.startswith("ID="):
                            gene_id = item.replace("ID=", "", 1)
                            break
                    if gene_id is None:
                        gene_id = gene
                else:
                    gene_id = gene

                return gene_id, f"{chrom}:{start}-{end}"

    raise ValueError(f"Gene not found in GFF: {gene}")


def get_parallel_workers(config: dict) -> int:
    """
    Return the number of parallel workers requested in the config.

    parallel:
      workers: 8

    Use 1 for sequential execution.
    """
    parallel_config = config.get("parallel", {})
    workers = parallel_config.get("workers", 1)

    if isinstance(workers, str):
        if workers.lower() == "auto":
            import os
            return max(1, os.cpu_count() or 1)
        workers = int(workers)

    workers = int(workers)

    if workers < 1:
        workers = 1

    return workers

def check_tool_runs(label: str, command: list) -> None:
    """
    Check that an external tool can actually run.

    This catches problems where the executable exists but fails at runtime,
    for example missing shared libraries.
    """
    import subprocess

    printable = " ".join(str(x) for x in command)

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"Tool '{label}' was found but failed to run.\n"
            f"Command: {printable}\n"
            f"Exit code: {result.returncode}\n"
            f"Message: {message}"
        )

    first_line = ""
    if result.stdout.strip():
        first_line = result.stdout.strip().splitlines()[0]
    elif result.stderr.strip():
        first_line = result.stderr.strip().splitlines()[0]

    if first_line:
        print(f"[OK] {label} runs: {first_line}")
    else:
        print(f"[OK] {label} runs")


def run_checks(config: dict) -> None:
    """
    Validate the PinPis configuration, input files, tools, and legacy scripts.

    This should run before any workflow step.
    """
    import shutil

    print("[pinpis] Checking configuration and inputs...")

    # ------------------------------------------------------------------
    # Project folders
    # ------------------------------------------------------------------
    output_dir = Path(require(config, "project.output_dir"))
    work_dir = Path(require(config, "project.work_dir"))
    tmp_dir = Path(require(config, "project.tmp_dir"))

    print(f"[OK] output_dir: {output_dir}")
    print(f"[OK] work_dir: {work_dir}")
    print(f"[OK] tmp_dir: {tmp_dir}")

    # ------------------------------------------------------------------
    # Required input files
    # ------------------------------------------------------------------
    input_paths = {
        "VCF": require(config, "inputs.vcf"),
        "GFF": require(config, "inputs.gff"),
        "reference genome": require(config, "inputs.reference_genome"),
        "samples ID file": require(config, "inputs.samples_id_file"),
        "sample years file": require(config, "inputs.sample_years_file"),
        "gene list": require(config, "inputs.gene_list"),
    }

    for label, path in input_paths.items():
        path_obj = Path(path)
        if not path_obj.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
        if path_obj.is_file() and path_obj.stat().st_size == 0:
            raise ValueError(f"{label} exists but is empty: {path}")
        print(f"[OK] {label}: {path}")

    # ------------------------------------------------------------------
    # Gene list
    # ------------------------------------------------------------------
    gene_list = Path(require(config, "inputs.gene_list"))
    genes = [line.strip() for line in gene_list.read_text().splitlines() if line.strip()]

    if len(genes) == 0:
        raise ValueError(f"Gene list is empty: {gene_list}")

    print(f"[OK] gene list contains {len(genes)} genes")

    # ------------------------------------------------------------------
    # Sampling groups
    # ------------------------------------------------------------------
    groups = require(config, "sampling.groups")

    if not isinstance(groups, list) or len(groups) < 2:
        raise ValueError("sampling.groups must contain at least two groups.")

    group_names = []

    for group in groups:
        name = group.get("name")
        population_file = group.get("population_file")
        output_dir_group = group.get("output_dir")

        if not name:
            raise ValueError("Each sampling group must have a name.")
        if not population_file:
            raise ValueError(f"Group {name} is missing population_file.")
        if not output_dir_group:
            raise ValueError(f"Group {name} is missing output_dir.")

        if name in group_names:
            raise ValueError(f"Duplicate group name: {name}")

        group_names.append(name)

        population_path = Path(population_file)
        if not population_path.exists():
            raise FileNotFoundError(
                f"Population file for group {name} not found: {population_file}"
            )
        if population_path.stat().st_size == 0:
            raise ValueError(
                f"Population file for group {name} is empty: {population_file}"
            )

        n_populations = sum(
            1 for line in population_path.read_text().splitlines()
            if line.strip()
        )

        print(
            f"[OK] group {name}: {population_file} "
            f"({n_populations} populations)"
        )

    # ------------------------------------------------------------------
    # External tools
    # ------------------------------------------------------------------
    external = config.get("external", {})


    required_tools = {
        "python": external.get("python", "python"),
        "bcftools": external.get("bcftools", "bcftools"),
        "samtools": external.get("samtools", "samtools"),
        "tabix": external.get("tabix", "tabix"),
        "bgzip": external.get("bgzip", "bgzip"),
        "yn00": external.get("yn00", "yn00"),
    }

    resolved_tools = {}

    for label, executable in required_tools.items():
        found = shutil.which(executable)

        # Allow absolute paths too.
        if found is None and Path(executable).exists():
            found = str(Path(executable))

        if found is None:
            raise FileNotFoundError(
                f"Required external tool not found for {label}: {executable}"
            )

        resolved_tools[label] = found
        print(f"[OK] tool {label}: {found}")

    # Actually run the tools. This catches broken binaries / missing libraries.
    check_tool_runs("python", [resolved_tools["python"], "--version"])
    check_tool_runs("bcftools", [resolved_tools["bcftools"], "--version"])
    check_tool_runs("samtools", [resolved_tools["samtools"], "--version"])
    check_tool_runs("tabix", [resolved_tools["tabix"], "--version"])
    check_tool_runs("bgzip", [resolved_tools["bgzip"], "--version"])

    # yn00 does not have a clean --version flag in all installations.
    # Running it with no arguments usually prints the PAML/YN00 banner and exits.

    yn00_result = subprocess.run(
        [resolved_tools["yn00"]],
        text=True,
        capture_output=True,
    )

    yn00_output = (yn00_result.stdout + "\n" + yn00_result.stderr).strip()

    if "ctl file open error" in yn00_output or "YN00" in yn00_output or "paml" in yn00_output.lower():
        print("[OK] yn00 runs: binary responds and expects a control file")
    else:
        raise RuntimeError(
            "Tool 'yn00' was found but did not respond as expected.\n"
            f"Command: {resolved_tools['yn00']}\n"
            f"Exit code: {yn00_result.returncode}\n"
            f"Message: {yn00_output}"
        )

    
    # GNU parallel is optional for the Python wrapper, but useful for legacy VCF cutting.
    parallel_exe = external.get("parallel", "parallel")
    parallel_found = shutil.which(parallel_exe)

    if parallel_found is None and Path(parallel_exe).exists():
        parallel_found = str(Path(parallel_exe))

    if parallel_found is None:
        print(
            f"[WARN] GNU parallel not found as '{parallel_exe}'. "
            "Some legacy steps may be slower."
        )
    else:
        print(f"[OK] tool parallel: {parallel_found}")
        check_tool_runs("parallel", [parallel_found, "--version"])    
        
        

    # ------------------------------------------------------------------
    # Legacy scripts
    # ------------------------------------------------------------------
    scripts_dir = Path(require(config, "legacy.scripts_dir"))

    if not scripts_dir.exists():
        raise FileNotFoundError(f"Legacy scripts directory not found: {scripts_dir}")

    required_scripts = [
        "multysampler.py",
        "sampler.py",
        "vcf_cutter_parallel.sh",
        "vcf_by_sample.sh",
        "vcf_to_cds.sh",
        "vcf_to_dna.sh",
        "CDS.py",
        "remove_indels.py",
        "dnds2.sh",
        "fasta_to_phylip_sequential.py",
        "bed_extractor.sh",
    ]

    for script in required_scripts:
        script_path = scripts_dir / script
        if not script_path.exists():
            raise FileNotFoundError(f"Required legacy script not found: {script_path}")
        if script_path.stat().st_size == 0:
            raise ValueError(f"Legacy script exists but is empty: {script_path}")

        print(f"[OK] legacy script: {script_path}")

    # ------------------------------------------------------------------
    # Parallel settings
    # ------------------------------------------------------------------
    workers = get_parallel_workers(config)
    print(f"[OK] parallel workers: {workers}")

    print("[pinpis] Check passed.")


def get_tool_versions(config: dict) -> dict:
    """
    Return short version strings for external tools.

    This is intentionally forgiving: if a tool version cannot be read,
    the summary records the error instead of crashing the pipeline.
    """
    import subprocess

    external = config.get("external", {})

    commands = {
        "python": [external.get("python", "python"), "--version"],
        "bcftools": [external.get("bcftools", "bcftools"), "--version"],
        "samtools": [external.get("samtools", "samtools"), "--version"],
        "tabix": [external.get("tabix", "tabix"), "--version"],
        "bgzip": [external.get("bgzip", "bgzip"), "--version"],
        "parallel": [external.get("parallel", "parallel"), "--version"],
    }

    versions = {}

    for name, command in commands.items():
        try:
            result = subprocess.run(
                command,
                text=True,
                capture_output=True,
                timeout=10,
            )

            text = (result.stdout or result.stderr or "").strip()
            first_line = text.splitlines()[0] if text else "unknown"

            if result.returncode == 0:
                versions[name] = first_line
            else:
                versions[name] = f"failed: {first_line}"

        except Exception as error:
            versions[name] = f"failed: {error}"

    # yn00 exits nonzero without a control file, but that still tells us the binary responds.
    yn00 = external.get("yn00", "yn00")
    try:
        result = subprocess.run(
            [yn00],
            text=True,
            capture_output=True,
            timeout=10,
        )

        text = (result.stdout + "\n" + result.stderr).strip()
        first_line = text.splitlines()[0] if text else "binary responds"
        versions["yn00"] = first_line

    except Exception as error:
        versions["yn00"] = f"failed: {error}"

    return versions


def write_run_summary(config: dict, command_name: str, config_path: str = None) -> None:
    """
    Write a human-readable reproducibility summary to results/run_summary.txt.
    """
    import datetime

    output_dir = Path(require(config, "project.output_dir"))
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_file = output_dir / "run_summary.txt"

    gene_list = Path(require(config, "inputs.gene_list"))
    genes = [
        line.strip()
        for line in gene_list.read_text().splitlines()
        if line.strip()
    ]

    groups = require(config, "sampling.groups")
    tool_versions = get_tool_versions(config)

    lines = []

    lines.append("PinPis run summary")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"Timestamp: {datetime.datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Project: {config.get('project', {}).get('name', 'unknown')}")
    lines.append(f"Command: {command_name}")
    lines.append(f"Config file: {config_path}")
    lines.append("")

    lines.append("Inputs")
    lines.append("-" * 80)
    lines.append(f"VCF: {require(config, 'inputs.vcf')}")
    lines.append(f"GFF: {require(config, 'inputs.gff')}")
    lines.append(f"Reference genome: {require(config, 'inputs.reference_genome')}")
    lines.append(f"Samples ID file: {require(config, 'inputs.samples_id_file')}")
    lines.append(f"Sample years file: {require(config, 'inputs.sample_years_file')}")
    lines.append(f"Gene list: {gene_list}")
    lines.append(f"Number of genes: {len(genes)}")
    lines.append("First genes:")
    for gene in genes[:10]:
        lines.append(f"  - {gene}")
    lines.append("")

    lines.append("Groups")
    lines.append("-" * 80)
    for group in groups:
        lines.append(f"Group: {group['name']}")
        lines.append(f"  Population file: {group['population_file']}")
        lines.append(f"  Output directory: {group['output_dir']}")
    lines.append("")

    lines.append("Parallelization")
    lines.append("-" * 80)
    lines.append(f"Workers: {get_parallel_workers(config)}")
    lines.append("")

    lines.append("Steps")
    lines.append("-" * 80)
    steps = config.get("steps", {})
    if steps:
        for step, enabled in steps.items():
            lines.append(f"{step}: {enabled}")
    else:
        lines.append("No explicit steps block found; default workflow behavior used.")
    lines.append("")

    lines.append("External tools")
    lines.append("-" * 80)
    external = config.get("external", {})
    for tool, command in external.items():
        lines.append(f"{tool}: {command}")
    lines.append("")

    lines.append("Tool versions")
    lines.append("-" * 80)
    for tool, version in tool_versions.items():
        lines.append(f"{tool}: {version}")
    lines.append("")

    lines.append("Outputs")
    lines.append("-" * 80)
    lines.append(f"Work directory: {require(config, 'project.work_dir')}")
    lines.append(f"Output directory: {require(config, 'project.output_dir')}")
    lines.append(f"Variant-site table: {output_dir / 'qc' / 'variant_sites_by_country.tsv'}")
    lines.append(f"Final table: {output_dir / 'pinpis_results.tsv'}")
    lines.append("")

    with summary_file.open("w") as handle:
        handle.write("\n".join(lines) + "\n")

    print(f"[pinpis] Run summary written to: {summary_file}")


def run_parallel_by_gene(config: dict, genes: list, task_function, task_name: str) -> None:
    """
    Run a per-gene task either sequentially or in parallel.

    task_function must accept:
        task_function(config, gene)
    """
    workers = get_parallel_workers(config)

    if workers == 1:
        print(f"[pinpis] Running {task_name} sequentially.")
        for gene in genes:
            task_function(config, gene)
        return

    print(f"[pinpis] Running {task_name} with {workers} parallel workers.")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_gene = {
            executor.submit(task_function, config, gene): gene
            for gene in genes
        }

        for future in as_completed(future_to_gene):
            gene = future_to_gene[future]

            try:
                future.result()
            except Exception as error:
                raise RuntimeError(
                    f"{task_name} failed for gene {gene}: {error}"
                ) from error

def run_cds_extraction_for_gene(config: dict, gene: str) -> None:
    scripts_dir = Path(require(config, "legacy.scripts_dir"))
    bed_extractor = scripts_dir / "bed_extractor.sh"
    vcf_to_cds = scripts_dir / "vcf_to_cds.sh"

    external = config.get("external", {})
    env = {
        "REFERENCE_GENOME": require(config, "inputs.reference_genome"),
        "BCFTOOLS_BIN": external.get("bcftools", "bcftools"),
        "SAMTOOLS_BIN": external.get("samtools", "samtools"),
        "BGZIP_BIN": external.get("bgzip", "bgzip"),
        "TABIX_BIN": external.get("tabix", "tabix"),
        "PYTHON_BIN": external.get("python", "python3"),
    }

    use_gene_name = str(config.get("annotation", {}).get("use_gene_name", False)).lower() in {
        "true", "1", "yes"
    }

    gff_file = require(config, "inputs.gff")

    genome_genes_dir = Path(
        config.get("annotation", {}).get("genome_genes_dir", "Genome/Genes")
    )
    genome_genes_dir.mkdir(parents=True, exist_ok=True)

    group_names = [group["name"] for group in require(config, "sampling.groups")]

    if len(group_names) < 2:
        raise ValueError("At least two groups are required for CDS extraction.")

    print(f"[pinpis] Extracting CDS for gene: {gene}")

    gene_id, genomic_region = get_gene_info_from_gff(
        gene=gene,
        gff_file=gff_file,
        use_gene_name=use_gene_name,
    )

    bed_file = genome_genes_dir / f"{gene_id}.bed"

    run_command([
        "bash",
        str(bed_extractor),
        gene_id,
        gff_file,
        str(genome_genes_dir),
    ], env=env)

    for group_name in group_names[:2]:
        group_dir = Path(require(config, "project.work_dir")) / gene / group_name

        if not group_dir.exists():
            raise FileNotFoundError(f"Expected group VCF directory not found: {group_dir}")

        run_command([
            "bash",
            str(vcf_to_cds),
            str(group_dir),
            genomic_region,
            str(bed_file),
        ], env=env)


def run_cds_extraction(config: dict) -> None:
    genes = read_gene_ids(config)
    run_parallel_by_gene(
        config=config,
        genes=genes,
        task_function=run_cds_extraction_for_gene,
        task_name="CDS extraction",
    )



def run_dnds_estimation_for_gene(config: dict, gene: str) -> None:
    scripts_dir = Path(require(config, "legacy.scripts_dir"))
    dnds_script = scripts_dir / "dnds2.sh"

    external = config.get("external", {})
    
    include_between = config.get("dnds", {}).get("include_between_groups", False)

    env = {
        "PYTHON_BIN": external.get("python", "python3"),
        "YN00_BIN": external.get("yn00", "yn00"),
        "INCLUDE_BETWEEN_GROUPS": str(include_between).lower(),
    }
    
    group_names = [group["name"] for group in require(config, "sampling.groups")]
    if len(group_names) < 2:
        raise ValueError("At least two groups are required for dN/dS estimation.")

    output_dir = Path(require(config, "project.output_dir"))
    output_dir.mkdir(parents=True, exist_ok=True)

    work_dir = Path(require(config, "project.work_dir"))

    print(f"[pinpis] Estimating dN/dS for gene: {gene}")

    group_1_dir = work_dir / gene / group_names[0]
    group_2_dir = work_dir / gene / group_names[1]

    if not group_1_dir.exists():
        raise FileNotFoundError(f"Missing group directory: {group_1_dir}")
    if not group_2_dir.exists():
        raise FileNotFoundError(f"Missing group directory: {group_2_dir}")

    output_file = output_dir / f"{gene}_dnds.txt"

    if output_file.exists():
        print(f"[pinpis] Removing existing dN/dS output: {output_file}")
        output_file.unlink()

    run_command([
        "bash",
        str(dnds_script),
        str(group_1_dir),
        str(group_2_dir),
        str(output_file),
        gene,
    ], env=env)


def run_dnds_estimation(config: dict) -> None:
    genes = read_gene_ids(config)
    run_parallel_by_gene(
        config=config,
        genes=genes,
        task_function=run_dnds_estimation_for_gene,
        task_name="dN/dS estimation",
    )
    

def read_fasta_sequences(fasta_path: Path) -> list:
    """
    Read sequences from a FASTA file into a list of uppercase strings.
    """
    sequences = []
    current = []

    with fasta_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)

        if current:
            sequences.append("".join(current).upper())

    return sequences


def count_variable_sites_in_sequences(sequences: list, label: str = "") -> int:
    """
    Count variable nucleotide sites across a list of aligned sequences.

    Only A/C/G/T are considered. Other symbols are ignored.
    """
    if not sequences:
        return 0

    length = len(sequences[0])

    for seq in sequences:
        if len(seq) != length:
            raise ValueError(f"Sequences have unequal lengths in {label}")

    variable_sites = 0

    for i in range(length):
        bases = set()

        for seq in sequences:
            base = seq[i]

            if base in "ACGT":
                bases.add(base)

                if len(bases) > 1:
                    variable_sites += 1
                    break

    return variable_sites


def count_variable_sites_in_fasta(fasta_path: Path) -> int:
    """
    Count variable sites in one aligned FASTA file.
    """
    sequences = read_fasta_sequences(fasta_path)
    return count_variable_sites_in_sequences(sequences, label=str(fasta_path))

def count_variable_sites_in_fasta(fasta_path: Path) -> int:
    """
    Count variable nucleotide sites in a FASTA alignment.

    Only A/C/G/T are considered. Other symbols are ignored for the
    purpose of deciding whether a column is variable.
    """
    from Bio import AlignIO

    alignment = AlignIO.read(str(fasta_path), "fasta")
    variable_sites = 0

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i].upper()
        bases = [base for base in column if base in "ACGT"]

        if len(set(bases)) > 1:
            variable_sites += 1

    return variable_sites


def count_sequences_in_fasta(fasta_path: Path) -> int:
    n_sequences = 0

    with fasta_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                n_sequences += 1

    return n_sequences



def collect_variant_site_rows_for_gene(config: dict, gene: str) -> list:
    work_dir = Path(require(config, "project.work_dir"))
    group_names = [group["name"] for group in require(config, "sampling.groups")]

    rows = []

    for group_name in group_names:
        seq_dir = work_dir / gene / group_name / "sequences"

        if not seq_dir.exists():
            print(f"[pinpis] Warning: missing sequence directory: {seq_dir}")
            continue

        for fasta_path in sorted(seq_dir.glob("*_CDS.fasta")):
            country = fasta_path.name.replace("_CDS.fasta", "")

            try:
                n_sequences = count_sequences_in_fasta(fasta_path)
                variant_sites = count_variable_sites_in_fasta(fasta_path)
            except Exception as error:
                print(f"[pinpis] Warning: failed to process {fasta_path}: {error}")
                n_sequences = "NA"
                variant_sites = "NA"

            rows.append({
                "gene": gene,
                "group": group_name,
                "country": country,
                "variant_sites": variant_sites,
                "n_sequences": n_sequences,
                "cds_file": str(fasta_path),
            })

    return rows


def collect_population_variant_sites_for_gene(config: dict, gene: str) -> list:
    """
    Count variant sites at the configured population/group level.

    Example:
        gene PF3D7_x, group Africa:
            combine all Africa country CDS FASTAs
            count variable sites across the combined alignment
    """
    work_dir = Path(require(config, "project.work_dir"))
    group_names = [group["name"] for group in require(config, "sampling.groups")]

    rows = []

    for group_name in group_names:
        seq_dir = work_dir / gene / group_name / "sequences"

        if not seq_dir.exists():
            print(f"[pinpis] Warning: missing sequence directory: {seq_dir}")
            continue

        fasta_files = sorted(seq_dir.glob("*_CDS.fasta"))

        if not fasta_files:
            print(f"[pinpis] Warning: no CDS FASTA files found in {seq_dir}")
            continue

        all_sequences = []
        country_files = []

        for fasta_path in fasta_files:
            sequences = read_fasta_sequences(fasta_path)
            all_sequences.extend(sequences)
            country_files.append(str(fasta_path))

        try:
            variant_sites = count_variable_sites_in_sequences(
                all_sequences,
                label=f"{gene}:{group_name}",
            )
        except Exception as error:
            print(
                f"[pinpis] Warning: failed population-level variant-site count "
                f"for {gene} {group_name}: {error}"
            )
            variant_sites = "NA"

        rows.append({
            "gene": gene,
            "population": group_name,
            "variant_sites": variant_sites,
            "n_sequences": len(all_sequences),
            "n_country_files": len(fasta_files),
        })

    return rows


def collect_country_variant_site_rows_for_gene(config: dict, gene: str) -> list:
    """
    Optional country/subpopulation-level variant-site QC.
    """
    work_dir = Path(require(config, "project.work_dir"))
    group_names = [group["name"] for group in require(config, "sampling.groups")]

    rows = []

    for group_name in group_names:
        seq_dir = work_dir / gene / group_name / "sequences"

        if not seq_dir.exists():
            print(f"[pinpis] Warning: missing sequence directory: {seq_dir}")
            continue

        for fasta_path in sorted(seq_dir.glob("*_CDS.fasta")):
            country = fasta_path.name.replace("_CDS.fasta", "")

            try:
                sequences = read_fasta_sequences(fasta_path)
                n_sequences = len(sequences)
                variant_sites = count_variable_sites_in_sequences(
                    sequences,
                    label=str(fasta_path),
                )
            except Exception as error:
                print(f"[pinpis] Warning: failed to process {fasta_path}: {error}")
                n_sequences = "NA"
                variant_sites = "NA"

            rows.append({
                "gene": gene,
                "group": group_name,
                "country": country,
                "variant_sites": variant_sites,
                "n_sequences": n_sequences,
                "cds_file": str(fasta_path),
            })

    return rows


def run_variant_site_counting(config: dict) -> None:
    """
    Count variant sites.

    Default:
        population-level counts only.

    Optional:
        country/subpopulation-level QC if variant_sites.write_country_qc is true.
    """
    output_dir = Path(require(config, "project.output_dir"))
    qc_dir = output_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    population_output = qc_dir / "variant_sites_by_population.tsv"
    country_output = qc_dir / "variant_sites_by_country.tsv"

    genes = read_gene_ids(config)
    workers = get_parallel_workers(config)

    write_country_qc = bool(
        config.get("variant_sites", {}).get("write_country_qc", False)
    )

    print(f"[pinpis] Writing population-level variant-site counts to: {population_output}")

    population_rows = []

    if workers == 1:
        print("[pinpis] Running population-level variant-site counting sequentially.")
        for gene in genes:
            population_rows.extend(collect_population_variant_sites_for_gene(config, gene))
    else:
        print(
            f"[pinpis] Running population-level variant-site counting "
            f"with {workers} parallel workers."
        )

        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_gene = {
                executor.submit(collect_population_variant_sites_for_gene, config, gene): gene
                for gene in genes
            }

            for future in as_completed(future_to_gene):
                gene = future_to_gene[future]
                try:
                    population_rows.extend(future.result())
                except Exception as error:
                    raise RuntimeError(
                        f"Population-level variant-site counting failed for gene {gene}: {error}"
                    ) from error

    population_rows = sorted(
        population_rows,
        key=lambda row: (row["gene"], row["population"]),
    )

    with population_output.open("w") as out:
        out.write("gene\tpopulation\tvariant_sites\tn_sequences\tn_country_files\n")

        for row in population_rows:
            out.write(
                f"{row['gene']}\t{row['population']}\t"
                f"{row['variant_sites']}\t{row['n_sequences']}\t"
                f"{row['n_country_files']}\n"
            )

    if write_country_qc:
        print(f"[pinpis] Writing country-level variant-site QC to: {country_output}")

        country_rows = []

        if workers == 1:
            print("[pinpis] Running country-level variant-site QC sequentially.")
            for gene in genes:
                country_rows.extend(collect_country_variant_site_rows_for_gene(config, gene))
        else:
            print(
                f"[pinpis] Running country-level variant-site QC "
                f"with {workers} parallel workers."
            )

            with ThreadPoolExecutor(max_workers=workers) as executor:
                future_to_gene = {
                    executor.submit(collect_country_variant_site_rows_for_gene, config, gene): gene
                    for gene in genes
                }

                for future in as_completed(future_to_gene):
                    gene = future_to_gene[future]
                    try:
                        country_rows.extend(future.result())
                    except Exception as error:
                        raise RuntimeError(
                            f"Country-level variant-site QC failed for gene {gene}: {error}"
                        ) from error

        country_rows = sorted(
            country_rows,
            key=lambda row: (row["gene"], row["group"], row["country"]),
        )

        with country_output.open("w") as out:
            out.write("gene\tgroup\tcountry\tvariant_sites\tn_sequences\tcds_file\n")

            for row in country_rows:
                out.write(
                    f"{row['gene']}\t{row['group']}\t{row['country']}\t"
                    f"{row['variant_sites']}\t{row['n_sequences']}\t"
                    f"{row['cds_file']}\n"
                )
    else:
        if country_output.exists():
            print(f"[pinpis] Removing old country-level QC table: {country_output}")
            country_output.unlink()

    print("[pinpis] Variant-site counting finished.")


def run_final_table(config: dict) -> None:
    """
    Merge dN/dS results with population-level variant-site counts and write
    the main PinPis result table.

    Output:
        results/pinpis_results.tsv
    """
    import pandas as pd

    output_dir = Path(require(config, "project.output_dir"))
    qc_dir = output_dir / "qc"

    variant_sites_file = qc_dir / "variant_sites_by_population.tsv"
    final_output = output_dir / "pinpis_results.tsv"

    if not variant_sites_file.exists():
        raise FileNotFoundError(
            f"Population-level variant-site table not found: {variant_sites_file}. "
            "Run count-sites first."
        )

    variant_df = pd.read_csv(variant_sites_file, sep="\t")

    required_variant_cols = {
        "gene", "population", "variant_sites", "n_sequences"
    }

    missing = required_variant_cols - set(variant_df.columns)
    if missing:
        raise ValueError(
            f"Population-level variant-site table is missing columns: {sorted(missing)}"
        )

    group_names = [group["name"] for group in require(config, "sampling.groups")]

    if len(group_names) < 2:
        raise ValueError("At least two groups are required for final table generation.")

    include_between_raw = config.get("dnds", {}).get("include_between_groups", False)

    if isinstance(include_between_raw, bool):
        include_between = include_between_raw
    else:
        include_between = str(include_between_raw).lower() in {"true", "1", "yes", "y"}


    all_rows = []

    for gene in read_gene_ids(config):
        dnds_file = output_dir / f"{gene}_dnds.txt"

        if not dnds_file.exists():
            print(f"[pinpis] Warning: missing dN/dS file for {gene}: {dnds_file}")
            continue

        dnds_df = pd.read_csv(dnds_file, sep="\t")

        expected_cols = ["population_1", "population_2", "dN", "dS", "dN_dS"]

        if not set(expected_cols).issubset(dnds_df.columns):
            dnds_df = pd.read_csv(
                dnds_file,
                sep="\t",
                header=None,
                names=expected_cols,
            )
            dnds_df = dnds_df[dnds_df["population_1"] != "population_1"]

        for col in ["dN", "dS", "dN_dS"]:
            dnds_df[col] = pd.to_numeric(dnds_df[col], errors="coerce")

        gene_variant_df = variant_df[variant_df["gene"] == gene].copy()

        # Country/subpopulation names are inferred from available CDS files.
        # We use them only to decide which dN/dS pairwise rows belong to each group.
        countries_by_group = {}

        for group in group_names:
            seq_dir = Path(require(config, "project.work_dir")) / gene / group / "sequences"
            countries = {
                fasta.name.replace("_CDS.fasta", "")
                for fasta in seq_dir.glob("*_CDS.fasta")
            }
            countries_by_group[group] = countries

        # Main within-population summaries.
        for group in group_names:
            countries = countries_by_group[group]

            group_dnds = dnds_df[
                dnds_df["population_1"].astype(str).isin(countries)
                & dnds_df["population_2"].astype(str).isin(countries)
            ].copy()

            group_variant = gene_variant_df[
                gene_variant_df["population"] == group
            ].copy()

            if group_dnds.empty:
                print(f"[pinpis] Warning: no dN/dS rows for {gene} {group}")
                continue

            if group_variant.empty:
                print(f"[pinpis] Warning: no variant-site row for {gene} {group}")
                variant_sites = "NA"
                n_sequences = "NA"
            else:
                variant_sites = group_variant["variant_sites"].iloc[0]
                n_sequences = group_variant["n_sequences"].iloc[0]

            all_rows.append({
                "gene": gene,
                "population": group,
                "dN": group_dnds["dN"].mean(),
                "dS": group_dnds["dS"].mean(),
                "dN_dS": group_dnds["dN_dS"].mean(),
                "variant_sites": variant_sites,
                "n_sequences": n_sequences,
                "n_pairwise_comparisons": len(group_dnds),
            })

        # Optional between-group summary.
        if include_between and len(group_names) >= 2:
            group_a = group_names[0]
            group_b = group_names[1]

            countries_a = countries_by_group[group_a]
            countries_b = countries_by_group[group_b]

            between_dnds = dnds_df[
                (
                    dnds_df["population_1"].astype(str).isin(countries_a)
                    & dnds_df["population_2"].astype(str).isin(countries_b)
                )
                |
                (
                    dnds_df["population_1"].astype(str).isin(countries_b)
                    & dnds_df["population_2"].astype(str).isin(countries_a)
                )
            ].copy()


            if not between_dnds.empty:
                all_rows.append({
                    "gene": gene,
                    "population": f"{group_a}_vs_{group_b}",
                    "dN": between_dnds["dN"].mean(),
                    "dS": between_dnds["dS"].mean(),
                    "dN_dS": between_dnds["dN_dS"].mean(),
                    "variant_sites": "NA",
                    "n_sequences": "NA",
                    "n_pairwise_comparisons": len(between_dnds),
                })

    final_df = pd.DataFrame(all_rows)

    if final_df.empty:
        raise ValueError("No rows were generated for the final table.")

    final_df.to_csv(final_output, sep="\t", index=False)

    print(f"[pinpis] Final population-level table written to: {final_output}")

        
def run_vcf_cutting(config: dict) -> None:
    scripts_dir = Path(require(config, "legacy.scripts_dir"))
    cutter = scripts_dir / "vcf_cutter_parallel.sh"

    use_gene_name = str(config.get("annotation", {}).get("use_gene_name", False))
    use_gene_name = "True" if use_gene_name.lower() in {"true", "1", "yes"} else "False"

    run_command([
        "bash",
        str(cutter),
        require(config, "inputs.gene_list"),
        require(config, "sampling.samples_file_dir"),
        require(config, "inputs.vcf"),
        require(config, "inputs.gff"),
        use_gene_name,
    ])


def step_enabled(config: dict, step_name: str, default: bool = True) -> bool:
    steps = config.get("steps", {})
    value = steps.get(step_name, default)

    if isinstance(value, bool):
        return value

    if isinstance(value, str):
        return value.lower() in {"true", "1", "yes", "y"}

    return bool(value)


def prepare_workflow(config: dict) -> None:
    """
    Shared setup used by all CLI subcommands.
    """
    run_checks(config)
    setup_directories(config)

def run_legacy_pipeline(config: dict) -> None:
    """
    Run the complete PinPis workflow.

    By default, all steps run. If a config contains a `steps:` block,
    individual steps can be disabled for development or reruns.
    """
    prepare_workflow(config)

    if step_enabled(config, "sampling", True):
        run_sampling(config)
    else:
        print("[pinpis] Skipping sampling step.")

    if step_enabled(config, "vcf_cutting", True):
        run_vcf_cutting(config)
    else:
        print("[pinpis] Skipping VCF cutting step.")

    if step_enabled(config, "cds_extraction", True):
        run_cds_extraction(config)
    else:
        print("[pinpis] Skipping CDS extraction step.")

    if step_enabled(config, "dnds_estimation", True):
        run_dnds_estimation(config)
    else:
        print("[pinpis] Skipping dN/dS estimation step.")

    if step_enabled(config, "variant_site_counting", True):
        run_variant_site_counting(config)
    else:
        print("[pinpis] Skipping variant-site counting step.")

    if step_enabled(config, "final_table", True):
        run_final_table(config)
    else:
        print("[pinpis] Skipping final table step.")

    print("[pinpis] Full workflow finished.")
