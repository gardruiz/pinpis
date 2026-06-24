import argparse
import sys

from .config import load_config, validate_config, ConfigError
from .commands import CommandError
from .workflow import (
    run_checks,
    prepare_workflow,
    write_run_summary,
    run_sampling,
    run_vcf_cutting,
    run_cds_extraction,
    run_dnds_estimation,
    run_variant_site_counting,
    run_final_table,
    run_legacy_pipeline,
)


def add_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML config file.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pinpis",
        description=(
            "PinPis: coding-sequence reconstruction and piN/piS or dN/dS "
            "pipeline from VCF data."
        ),
    )
    
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    check_parser = subparsers.add_parser(
        "check",
        help="Validate config, input files, external tools, and legacy scripts.",
    )
    add_config_argument(check_parser)
    
    run_parser = subparsers.add_parser(
        "run",
        help="Run the workflow using the steps enabled in the YAML config.",
    )
    add_config_argument(run_parser)

    sample_parser = subparsers.add_parser(
        "sample",
        help="Select sample IDs for each configured population group.",
    )
    add_config_argument(sample_parser)

    cut_vcf_parser = subparsers.add_parser(
        "cut-vcf",
        help="Create gene-level and population-level VCF files.",
    )
    add_config_argument(cut_vcf_parser)

    extract_cds_parser = subparsers.add_parser(
        "extract-cds",
        help="Generate per-sample CDS FASTA files from VCFs.",
    )
    add_config_argument(extract_cds_parser)

    estimate_parser = subparsers.add_parser(
        "estimate",
        help="Run dN/dS estimation with yn00.",
    )
    add_config_argument(estimate_parser)

    count_sites_parser = subparsers.add_parser(
        "count-sites",
        help="Count variable sites in each country/population CDS FASTA file.",
    )
    add_config_argument(count_sites_parser)

    merge_parser = subparsers.add_parser(
        "merge",
        help="Merge dN/dS and variant-site tables into final PinPis output.",
    )
    add_config_argument(merge_parser)

    return parser


def load_and_validate(config_path: str) -> dict:
    config = load_config(config_path)
    validate_config(config)
    return config


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        config = load_and_validate(args.config)

        if args.command == "check":
            run_checks(config)

        elif args.command == "run":
            run_legacy_pipeline(config)

        elif args.command == "sample":
            prepare_workflow(config)
            run_sampling(config)

        elif args.command == "cut-vcf":
            prepare_workflow(config)
            run_vcf_cutting(config)

        elif args.command == "extract-cds":
            prepare_workflow(config)
            run_cds_extraction(config)

        elif args.command == "estimate":
            prepare_workflow(config)
            run_dnds_estimation(config)

        elif args.command == "count-sites":
            prepare_workflow(config)
            run_variant_site_counting(config)

        elif args.command == "merge":
            prepare_workflow(config)
            run_final_table(config)

        else:
            parser.error(f"Unknown command: {args.command}")
       
        if args.command != "check":
            write_run_summary(
                config=config,
                command_name=args.command,
                config_path=args.config,
            )

    except (ConfigError, FileNotFoundError, RuntimeError, CommandError, ValueError) as error:
        print(f"[pinpis] Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
