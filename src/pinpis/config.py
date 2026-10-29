from pathlib import Path
import yaml
from typing import Union


class ConfigError(Exception):
    pass


def load_config(config_path: Union[str, Path]) -> dict:
    config_path = Path(config_path)

    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")

    with config_path.open("r") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ConfigError("Config file is empty or invalid.")

    return config


def require(config: dict, path: str):
    """
    Get a nested config value using dot notation.

    Example:
        require(config, "inputs.vcf")
    """
    current = config
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            raise ConfigError(f"Missing required config value: {path}")
        current = current[key]
    return current


def validate_config(config: dict) -> None:
    required_paths = [
        "project.output_dir",
        "project.work_dir",
        "project.tmp_dir",
        "inputs.vcf",
        "inputs.gff",
        "inputs.reference_genome",
        "inputs.samples_id_file",
        "inputs.sample_years_file",
        "inputs.gene_list",
        "sampling.number_of_samples",
        "sampling.samples_file_dir",
        "sampling.groups",
        "legacy.scripts_dir",
    ]

    for path in required_paths:
        require(config, path)
    
    # Validate input files
    input_files = [
        "inputs.vcf",
        "inputs.gff",
        "inputs.reference_genome",
    ]

    for key in input_files:
        value = require(config, key)

        if not isinstance(value, str) or not value.strip():
            raise ConfigError(
                f"Missing or empty input file path: {key}"
            )

        file_path = Path(value).expanduser()

        if not file_path.is_file():
            raise ConfigError(
                f"Input file not found: {key} = {file_path}"
            )

    groups = require(config, "sampling.groups")
    if not isinstance(groups, list) or len(groups) < 2:
        raise ConfigError("sampling.groups must contain at least two groups.")

    for group in groups:
        if "name" not in group:
            raise ConfigError("Each sampling group requires a name.")
        if "population_file" not in group:
            raise ConfigError(f"Group {group.get('name', '<unknown>')} is missing population_file.")
        if "output_dir" not in group:
            raise ConfigError(f"Group {group.get('name', '<unknown>')} is missing output_dir.")
