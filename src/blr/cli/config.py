"""
Update configuration file. If no --set option is given the current settings are printed.
"""
import sys
import os
import logging
from ruamel.yaml import YAML
from snakemake.utils import validate
import pkg_resources

logger = logging.getLogger(__name__)
DEFAULT_PATH = "blr.yaml"
SCHEMA_FILE = "config.schema.yaml"


def main(args):
    # Script is based on repos NBISSweden/IgDisover config script.
    # Link https://github.com/NBISweden/IgDiscover/blob/master/src/igdiscover/cli/config.py

    if args.set:
        change_config(args.file, args.set)
    else:
        configs, yaml = load_yaml(args.file)
        print(f"--- CONFIGS IN: {args.file} ---")
        yaml.dump(configs, stream=sys.stdout)


def change_config(filename, changes_set):
    """
    Change config YAML file at filename using the changes_set key-value pairs.
    :param filename: string with path to YAML config file to change.
    :param changes_set: dict with changes to incorporate.
    """
    # Get configs from file.
    configs, yaml = load_yaml(filename)

    # Update configs
    for key, value in changes_set:
        if key in configs:
            value = YAML(typ='safe').load(value)
            logger.info(f"Changing value of '{key}': {configs[key]} --> {value}.")
            configs[key] = value
        else:
            logger.warning(f"KEY = {key} not in config. Config not updated with set ({key}, {value})")

    # Confirm that configs is valid.
    schema_path = pkg_resources.resource_filename("blr", SCHEMA_FILE)
    validate(configs, schema_path)

    # Write first to temporary file then overwrite filename.
    tmpfile = filename + ".tmp"
    with open(tmpfile, "w") as file:
        yaml.dump(configs, stream=file)
    os.rename(tmpfile, filename)


def load_yaml(filename):
    """
    Load YAML file and return the yaml object and data.
    :param filename: Path to YAML file
    :return: (data, yaml).
    """
    with open(filename) as file:
        yaml = YAML()
        data = yaml.load(file)
    return data, yaml


def add_arguments(parser):
    parser.add_argument("--set", nargs=2, metavar=("KEY", "VALUE"), action="append",
                        help="Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys. For empty values "
                             "write 'null'. Can be given multiple times.")
    parser.add_argument("--file", default=DEFAULT_PATH,
                        help="Configuration file to modify. Default: %(default)s in current directory.")
