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

    # Get configs from file.
    yaml_config = YAML()
    with open(args.file) as configfile:
        configs = yaml_config.load(configfile)

    # Udpate configs with SET or print current configs.
    if args.set:
        # Get schema form package data
        schema_path = pkg_resources.resource_filename("blr", SCHEMA_FILE)
        yaml_schema = YAML(typ="safe")
        with open(schema_path) as schemafile:
            schema = yaml_schema.load(schemafile)

        # Update configs
        for key, value in args.set:
            if key in schema['properties']:
                configs[key] = update_value(key, value, schema)
            else:
                logging.warning(f"KEY = {key} not in schema. Config not updated with set ({key}, {value})")

        # Confirm that data is valid.
        validate(configs, schema_path)

        # Write first to tmp then overwrite target.
        tmpfile = args.file + ".tmp"
        with open(tmpfile, "w") as file:
            yaml_config.dump(configs, stream=file)
        os.rename(tmpfile, args.file)
    else:
        print(f"--------Configs in: {args.file}--------")
        yaml_config.dump(configs, stream=sys.stdout)


def update_value(key, value, schema):
    """
    Update value to match allowed type in schema.
    :param key: Strint with name of parameter to update
    :param value: String with value to updates
    :param schema: dict with type information about parameter.
    :return: value
    """
    properties = schema['properties'][key]

    print(properties["type"])
    print(value == "null")
    print("null" in properties["type"])
    if value == "null" and "null" in properties["type"]:
        print("HEAF")
        return None

    if "integer" in properties["type"]:
        return int(value)

    if "number" in properties["type"]:
        return float(value)

    return value


def add_arguments(parser):
    arg = parser.add_argument
    arg("--set", nargs=2, default=[], metavar=("KEY", "VALUE"), action="append",
        help="Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys. To leave VALUE empty write 'null'"
             "Can be given multiple times.")
    arg("--file", default=DEFAULT_PATH,
        help="Configuration file to modify. Default: %(default)s in current directory.")
