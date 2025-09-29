import os
from collections.abc import MutableMapping as DictMixin
from importlib.resources import files
from shutil import copyfile

from appdirs import user_config_dir
from loguru import logger
from yaml import dump, full_load

__all__ = ["config", "manage_config"]


class Config(DictMixin):
    """Configuration class.

    State will be shared across config objects.
    """

    # Store shared state, Borg pattern
    __shared_state = {}
    config_dir = user_config_dir("genomepy")
    config_file = os.path.join(config_dir, "genomepy.yaml")
    default_config_file = str(files("genomepy.config").joinpath("default.yaml"))

    def __init__(self):
        """
        Create a Config object and read a config file.
        """
        self.__dict__ = self.__shared_state

        # Lookup the config file according to XDG hierarchy
        if not os.path.exists(self.config_dir):
            os.makedirs(self.config_dir)
        if not os.path.exists(self.config_file):
            copyfile(self.default_config_file, self.config_file)
        self.config = {}
        self.load(self.config_file)

    def load(self, path=None):
        """
        Load yaml-formatted config file.

        Parameters
        ----------
        path : str | None
            path to config file
        """
        if path is None:
            path = self.config_file
        with open(path) as f:
            self.config = full_load(f)
            if self.config is None:
                logger.warning("config file is empty!")
                self.config = {}

    def save(self, path=None):
        """
        Save current state of config dictionary.

        Parameters
        ----------
        path : str | None
            path to config file
        """
        if path is None:
            path = self.config_file
        with open(path, "w") as f:
            f.write(dump(self.config, default_flow_style=False))

    def __getitem__(self, key):
        return self.config.__getitem__(key)

    def __delitem__(self, key):
        self.config.__delitem__(key)

    def __setitem__(self, key, value):
        return self.config.__setitem__(key, value)

    def __len__(self):
        return self.config.__len__()

    def __iter__(self):
        return self.config.__iter__()

    def __str__(self):
        return (
            "{"
            + ", ".join(
                ": ".join([str(key), str(value)]) for key, value in self.items()
            )
            + "}"
        )

    def keys(self):
        return self.config.keys()


config = Config()


def generate_config():
    # overwrite existing config
    copyfile(config.default_config_file, config.config_file)
    config.load()
    logger.info(f"Created config file {config.config_file}")


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
