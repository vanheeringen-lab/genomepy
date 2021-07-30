import os
from shutil import copyfile

from appdirs import user_config_dir
from norns import config as cfg

__all__ = ["config", "manage_config"]

config = cfg("genomepy", default="config/default.yaml")


def generate_config():
    config_dir = user_config_dir("genomepy")
    new_config = os.path.join(config_dir, "genomepy.yaml")

    # existing config must be removed before norns picks up the default again
    os.makedirs(config_dir, exist_ok=True)
    if os.path.exists(new_config):
        os.remove(new_config)

    default_config = cfg("genomepy", default="config/default.yaml").config_file
    copyfile(default_config, new_config)
    config.config_file = new_config
    print(f"Created config file {new_config}")


def manage_config(command):
    """
    Manage the genomepy configuration

    Parameters
    ----------
    command : str
        command to perform. Options:

        file
            return config filepath
        show
            return config content
        generate
            create new config file
    """
    if command in ["file", "path"]:
        print(config.config_file)

    elif command in ["show", "list"]:
        with open(config.config_file) as f:
            print(f.read())

    elif command == "generate":
        generate_config()

    else:
        raise ValueError(
            f"Invalid config command: '{command}'. "
            "Options: 'file', 'show' or 'generate'."
        )
