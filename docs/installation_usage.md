# PinPis installation and usage guide

PinPis is an internal reproducible pipeline for estimating population-level coding-sequence dN/dS or Pin/Pis-like summaries from VCF data.

The current version was developed for *Plasmodium falciparum* Pf7 analyses. It reconstructs coding sequences from VCF data, estimates dN/dS using PAML `yn00`, counts population-level variant sites, and produces a clean gene-by-population result table.

This pipeline is intended for reproducible analysis of the *P. falciparum* project. It is not yet intended as a general public software package.

---

## 1. Installation

### 1.1 Enter the PinPis directory

```bash
cd /path/to/pinpis
```

The project directory should contain:

```text
pinpis/
├── environment.yml
├── pyproject.toml
├── README.md
├── config/
├── docs/
├── scripts/
│   └── legacy/
└── src/
    └── pinpis/
```

### 1.2 Create the conda environment

```bash
conda env create -f environment.yml
conda activate pinpis
```

The environment should include:

```text
bcftools 1.21
samtools
tabix
bgzip
PAML / yn00
GNU parallel
Python
Biopython
pandas
PyYAML
```

Check the installed tools:

```bash
bcftools --version
samtools --version
yn00
```

`yn00` may print an error about a missing control file when run without arguments. That is normal; it means the binary is available and expects a control file.

### 1.3 Install PinPis in editable mode

If the package was not automatically installed by the environment file, run:

```bash
python -m pip install -e .
```

Check that the command is available:

```bash
pinpis --help
```

---

## 2. Basic usage

PinPis is controlled by a YAML config file.

Before running the pipeline, always check the config:

```bash
pinpis check --config config/example_pf7.yaml
```

Run the full pipeline:

```bash
pinpis run --config config/example_pf7.yaml
```

Or run individual steps:

```bash
pinpis sample --config config/example_pf7.yaml
pinpis cut-vcf --config config/example_pf7.yaml
pinpis extract-cds --config config/example_pf7.yaml
pinpis estimate --config config/example_pf7.yaml
pinpis count-sites --config config/example_pf7.yaml
pinpis merge --config config/example_pf7.yaml
```

Recommended order for a fresh run:

```bash
pinpis check --config config/example_pf7.yaml
pinpis sample --config config/example_pf7.yaml
pinpis cut-vcf --config config/example_pf7.yaml
pinpis extract-cds --config config/example_pf7.yaml
pinpis estimate --config config/example_pf7.yaml
pinpis count-sites --config config/example_pf7.yaml
pinpis merge --config config/example_pf7.yaml
```

The final result table is:

```text
results/pinpis_results.tsv
```

---

## 3. Workflow overview

The pipeline performs these steps:

```text
sample IDs
→ cut VCF by gene and population group
→ reconstruct CDS FASTA files
→ estimate dN/dS with yn00
→ count population-level variant sites
→ merge final result table
```

### 3.1 `sample`

```bash
pinpis sample --config config/example_pf7.yaml
```

Selects sample IDs for each population group and subpopulation/country.

Example output:

```text
Samples/SamplesID/Africa/benin.txt
Samples/SamplesID/Africa/ghana.txt
Samples/SamplesID/Asia/cambodia.txt
```

### 3.2 `cut-vcf`

```bash
pinpis cut-vcf --config config/example_pf7.yaml
```

Creates gene-level VCFs and population/subpopulation-specific VCFs.

Example output:

```text
variants/<gene>/<gene>.vcf.gz
variants/<gene>/Africa/benin.vcf.gz
variants/<gene>/Asia/cambodia.vcf.gz
```

### 3.3 `extract-cds`

```bash
pinpis extract-cds --config config/example_pf7.yaml
```

Reconstructs coding sequences for each gene and population group.

Example output:

```text
variants/<gene>/Africa/sequences/benin_CDS.fasta
variants/<gene>/Asia/sequences/cambodia_CDS.fasta
```

### 3.4 `estimate`

```bash
pinpis estimate --config config/example_pf7.yaml
```

Runs PAML `yn00` to estimate dN, dS, and dN/dS.

Detailed gene-level output:

```text
results/<gene>_dnds.txt
results/<gene>/*.out
```

### 3.5 `count-sites`

```bash
pinpis count-sites --config config/example_pf7.yaml
```

Counts variant sites at the population level.

Default output:

```text
results/qc/variant_sites_by_population.tsv
```

Optional country-level QC output:

```text
results/qc/variant_sites_by_country.tsv
```

### 3.6 `merge`

```bash
pinpis merge --config config/example_pf7.yaml
```

Creates the main output table:

```text
results/pinpis_results.tsv
```

---

## 4. Main output

The main output is:

```text
results/pinpis_results.tsv
```

Columns:

| Column | Description |
|---|---|
| `gene` | Gene ID |
| `population` | Population/group name from the config |
| `dN` | Mean dN across relevant pairwise subpopulation comparisons |
| `dS` | Mean dS across relevant pairwise subpopulation comparisons |
| `dN_dS` | Mean dN/dS across relevant pairwise subpopulation comparisons |
| `variant_sites` | Population-level variant-site count |
| `n_sequences` | Number of CDS sequences used for the population-level variant-site count |
| `n_pairwise_comparisons` | Number of pairwise subpopulation comparisons summarized |

Example:

```text
gene            population  dN        dS        dN_dS     variant_sites  n_sequences  n_pairwise_comparisons
PF3D7_1303800   Africa      0.00088   0.00124   0.71498   1443           347          28
PF3D7_1303800   Asia        0.00053   0.00071   0.77731   481            350          28
```

If between-population comparisons are enabled, additional rows are added:

```text
PF3D7_1303800   Africa_vs_Asia   ...
```

For between-population rows, `variant_sites` and `n_sequences` are set to `NA`, because variant-site counts are population-specific.

---

## 5. The config YAML file

The config file controls the entire pipeline.

Example:

```yaml
project:
  name: pf7_africa_asia
  output_dir: results
  work_dir: variants
  tmp_dir: tmp

inputs:
  vcf: /path/to/Pf3D7_filtered.vcf.gz
  gff: /path/to/Pfalciparum.gff
  reference_genome: /path/to/Pfalciparum.genome.fasta
  samples_id_file: /path/to/Pf7_samples.txt
  sample_years_file: /path/to/Sample_years.txt
  gene_list: /path/to/gene_ids.txt

sampling:
  number_of_samples: 50
  samples_file_dir: Samples/SamplesID/
  groups:
    - name: Africa
      population_file: /path/to/country_list_africa.txt
      output_dir: Samples/SamplesID/Africa
    - name: Asia
      population_file: /path/to/country_list_asia.txt
      output_dir: Samples/SamplesID/Asia

annotation:
  use_gene_name: false

pipelines:
  - dnds

dnds:
  include_between_groups: false

variant_sites:
  write_country_qc: false

external:
  python: python
  bcftools: bcftools
  samtools: samtools
  tabix: tabix
  bgzip: bgzip
  yn00: yn00
  parallel: parallel

parallel:
  workers: 10

legacy:
  scripts_dir: scripts/legacy

steps:
  sampling: true
  vcf_cutting: true
  cds_extraction: true
  dnds_estimation: true
  variant_site_counting: true
  final_table: true
```

### 5.1 `project`

```yaml
project:
  name: pf7_africa_asia
  output_dir: results
  work_dir: variants
  tmp_dir: tmp
```

| Field | Description |
|---|---|
| `name` | Project name recorded in `results/run_summary.txt` |
| `output_dir` | Directory where result tables are written |
| `work_dir` | Directory where gene VCFs and CDS FASTA files are written |
| `tmp_dir` | Temporary directory |

Recommended values:

```yaml
output_dir: results
work_dir: variants
tmp_dir: tmp
```

### 5.2 `inputs`

```yaml
inputs:
  vcf: /path/to/Pf3D7_filtered.vcf.gz
  gff: /path/to/Pfalciparum.gff
  reference_genome: /path/to/Pfalciparum.genome.fasta
  samples_id_file: /path/to/Pf7_samples.txt
  sample_years_file: /path/to/Sample_years.txt
  gene_list: /path/to/gene_ids.txt
```

| Field | Description |
|---|---|
| `vcf` | Input VCF file. Must be bgzipped and indexed |
| `gff` | Genome annotation file |
| `reference_genome` | Reference FASTA used for normalization and consensus |
| `samples_id_file` | Sample metadata file used by the sampling scripts |
| `sample_years_file` | Sample year/QC file |
| `gene_list` | Text file with one gene ID per line |

The gene list should look like:

```text
PF3D7_1303800
PF3D7_0522400
PF3D7_0511500
```

### 5.3 `sampling`

```yaml
sampling:
  number_of_samples: 50
  samples_file_dir: Samples/SamplesID/
  groups:
    - name: Africa
      population_file: /path/to/country_list_africa.txt
      output_dir: Samples/SamplesID/Africa
    - name: Asia
      population_file: /path/to/country_list_asia.txt
      output_dir: Samples/SamplesID/Asia
```

| Field | Description |
|---|---|
| `number_of_samples` | Number of samples to select per subpopulation/country |
| `samples_file_dir` | Parent directory where sample ID files are stored |
| `groups` | List of population groups to compare |

Each group requires:

| Field | Description |
|---|---|
| `name` | Population/group name used in outputs |
| `population_file` | File listing countries/subpopulations in that group |
| `output_dir` | Directory where sample ID files for that group are written |

Example population file:

```text
benin
gambia
ghana
kenya
mali
senegal
sudan
```

For the current Pf7 analysis, the configured groups are continents:

```text
Africa
Asia
```

However, the pipeline uses the generic term `population`, because users may later define groups by region, country, treatment status, or other categories.

### 5.4 `annotation`

```yaml
annotation:
  use_gene_name: false
```

If `use_gene_name` is `false`, the pipeline searches the GFF using gene IDs such as:

```text
ID=PF3D7_1303800
```

If `use_gene_name` is `true`, it searches using gene names instead.

For the current Pf7 analysis, use:

```yaml
use_gene_name: false
```

### 5.5 `dnds`

```yaml
dnds:
  include_between_groups: false
```

Controls whether between-population comparisons are estimated.

If `false`, the pipeline estimates only within-population comparisons.

For 7 Africa countries and 7 Asia countries:

```text
within Africa: 28 comparisons
within Asia:   28 comparisons
total:         56 comparisons per gene
```

If `true`, the pipeline also estimates between-population comparisons:

```text
Africa vs Asia: 49 comparisons
total:          105 comparisons per gene
```

When `include_between_groups: true`, the final table includes rows such as:

```text
Africa_vs_Asia
```

Recommended default:

```yaml
dnds:
  include_between_groups: false
```

Use `true` only when cross-population pairwise estimates are needed.

### 5.6 `variant_sites`

```yaml
variant_sites:
  write_country_qc: false
```

Controls variant-site counting outputs.

Default behavior:

```yaml
write_country_qc: false
```

This writes only:

```text
results/qc/variant_sites_by_population.tsv
```

Population-level variant sites are calculated by combining all CDS sequences for a gene within a population group and counting alignment columns with more than one A/C/G/T base.

Optional behavior:

```yaml
write_country_qc: true
```

This also writes:

```text
results/qc/variant_sites_by_country.tsv
```

Country-level QC is useful for debugging, but it is slower and is not required for the main final table.

### 5.7 `external`

```yaml
external:
  python: python
  bcftools: bcftools
  samtools: samtools
  tabix: tabix
  bgzip: bgzip
  yn00: yn00
  parallel: parallel
```

These are the external commands used by the pipeline.

For reproducibility, do not hard-code machine-specific paths in the example config. Instead, activate the conda environment and use command names:

```yaml
bcftools: bcftools
samtools: samtools
yn00: yn00
```

`pinpis check` verifies that these tools exist and can run.

### 5.8 `parallel`

```yaml
parallel:
  workers: 10
```

Controls the number of gene-level workers.

Use fewer workers for debugging:

```yaml
parallel:
  workers: 1
```

Use more workers for larger runs:

```yaml
parallel:
  workers: 4
```

or:

```yaml
parallel:
  workers: 10
```

Do not set this too high before testing, because dN/dS estimation launches many external `yn00` jobs.

### 5.9 `legacy`

```yaml
legacy:
  scripts_dir: scripts/legacy
```

Directory containing the legacy scripts wrapped by the package.

This should normally remain:

```yaml
scripts_dir: scripts/legacy
```

### 5.10 `steps`

```yaml
steps:
  sampling: true
  vcf_cutting: true
  cds_extraction: true
  dnds_estimation: true
  variant_site_counting: true
  final_table: true
```

This controls what is run by:

```bash
pinpis run --config config/example_pf7.yaml
```

For example, to rerun only dN/dS estimation and merging:

```yaml
steps:
  sampling: false
  vcf_cutting: false
  cds_extraction: false
  dnds_estimation: true
  variant_site_counting: false
  final_table: true
```

For individual stages, it is usually clearer to use subcommands:

```bash
pinpis estimate --config config/example_pf7.yaml
pinpis merge --config config/example_pf7.yaml
```

---

## 6. Reproducibility outputs

Every successful command except `check` writes:

```text
results/run_summary.txt
```

This records:

```text
timestamp
command
config file
input files
gene list
number of genes
population groups
parallel workers
external tools
tool versions
expected outputs
```

Keep this file with the final results used for analysis.

---

## 7. Common checks

Check final table:

```bash
head results/pinpis_results.tsv
wc -l results/pinpis_results.tsv
```

Check population-level variant sites:

```bash
column -t results/qc/variant_sites_by_population.tsv | less -S
```

Check detailed dN/dS rows per gene:

```bash
gene="PF3D7_1303800"
wc -l results/${gene}_dnds.txt
```

For 7 subpopulations in each of two groups:

```text
57 lines  = within-population only
106 lines = within + between-population comparisons
```

Check between-population rows:

```bash
grep "Africa_vs_Asia" results/pinpis_results.tsv | head
```

---

## 8. Troubleshooting

### 8.1 `bcftools` library error

Example:

```text
bcftools: error while loading shared libraries: libcrypto.so.1.0.0
```

This means the `bcftools` executable exists but cannot run in the current environment.

Use the provided environment:

```bash
conda env create -f environment.yml
conda activate pinpis
pinpis check --config config/example_pf7.yaml
```

### 8.2 No CDS FASTA files found

Example:

```text
Error: no CDS FASTA files found in variants/<gene>/Africa
```

This means `estimate` was run before successful CDS extraction.

Run:

```bash
pinpis cut-vcf --config config/example_pf7.yaml
pinpis extract-cds --config config/example_pf7.yaml
```

Then check:

```bash
find variants/<gene> -name "*_CDS.fasta"
```

### 8.3 Between-population rows are missing

Make sure the config has:

```yaml
dnds:
  include_between_groups: true
```

Then rerun:

```bash
pinpis estimate --config config/example_pf7.yaml
pinpis merge --config config/example_pf7.yaml
```

### 8.4 Country-level variant-site QC is slow

Use:

```yaml
variant_sites:
  write_country_qc: false
```

This keeps only the faster population-level variant-site counting.
