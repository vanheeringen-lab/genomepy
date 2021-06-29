import os
from shutil import copyfile

from appdirs import user_config_dir
from loguru import logger
from norns import config as cfg

config = cfg("genomepy", default="cfg/default.yaml")


def generate_config():
    config_dir = user_config_dir("genomepy")
    new_config = os.path.join(config_dir, "genomepy.yaml")

    # existing config must be removed before norns picks up the default again
    os.makedirs(config_dir, exist_ok=True)
    if os.path.exists(new_config):
        os.remove(new_config)

    default_config = cfg("genomepy", default="cfg/default.yaml").config_file
    copyfile(default_config, new_config)
    config.config_file = new_config
    logger.info(f"Created config file {new_config}")


def manage_config(cmd):
    """Manage the genomepy config file."""
    if cmd in ["file", "path"]:
        print(config.config_file)

    elif cmd in ["show", "list"]:
        with open(config.config_file) as f:
            print(f.read())

    elif cmd == "generate":
        generate_config()

    else:
        logger.error(
            f"Invalid config command: '{cmd}'. Options: 'file', 'show' or 'generate'."
        )
        os._exit(0)  # noqa: error was caught, now exit silently.
