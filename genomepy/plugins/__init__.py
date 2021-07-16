"""Plugin class, modules & related functions"""
import os
import re

from genomepy.config import config

__all__ = ["Plugin", "manage_plugins", "get_active_plugins"]


class Plugin:
    """Plugin base class."""

    def __init__(self):
        self.name = convert(type(self).__name__).replace("_plugin", "")
        self.active = False

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False

    def after_genome_download(self, genome, threads, force):
        raise NotImplementedError("plugin should implement this method")

    def get_properties(self, genome):
        raise NotImplementedError("plugin should implement this method")


def convert(name: str) -> str:
    """
    Convert CamelCase to underscore. e.g. StarPlugin -> star_plugin

    Parameters
    ----------
    name : str
        Camelcase string

    Returns
    -------
    name : str
        Converted name
    """
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def list_plugins() -> list:
    plugin_dir = os.path.dirname(os.path.realpath(__file__))
    plugin_files = [f for f in os.listdir(plugin_dir) if f.endswith(".py")]
    plugin_names = [f[:-3] for f in plugin_files if not f.startswith("_")]
    return plugin_names


def init_plugins():
    """
    create a dictionary of plugin instances

    Returns
    -------
    plugins : dictionary
        key is plugin name, value Plugin object
    """
    # import plugins
    for plugin in list_plugins():
        __import__(f"genomepy.plugins.{plugin}")

    # for each Plugin subclass, save an instance to a dict
    d = {}
    active_plugins = config.get("plugin", [])
    for c in Plugin.__subclasses__():
        ins = c()

        if ins.name in active_plugins:
            ins.activate()

        d[ins.name] = ins

    return d


PLUGINS = init_plugins()


def get_active_plugins() -> list:
    """Returns all active plugin instances."""
    return [inst for name, inst in PLUGINS.items() if inst.active]


def activate(name):
    """Activate plugin.

    Parameters
    ----------
    name : str
        Plugin name.
    """
    if name in PLUGINS:
        PLUGINS[name].activate()
    else:
        raise ValueError(f"plugin {name} not found")


def deactivate(name):
    """Deactivate plugin.

    Parameters
    ----------
    name : str
        Plugin name.
    """
    if name in PLUGINS:
        PLUGINS[name].deactivate()
    else:
        raise ValueError(f"plugin {name} not found")


def show_plugins():
    active_plugins = config.get("plugin", [])
    print("{:20}{}".format("plugin", "enabled"))
    for plugin in sorted(PLUGINS):
        print(
            "{:20}{}".format(plugin, {False: "", True: "*"}[plugin in active_plugins])
        )


def manage_plugins(command: str, plugin_names: list = None):
    """
    Manage genomepy plugins

    Parameters
    ----------
    command : str
        command to perform. Options:

        list
            show plugins and status
        enable
            enable plugins
        disable
            disable plugins
    plugin_names : list
        plugin names for the enable/disable command
    """
    if command in ["show", "list"]:
        return show_plugins()

    active_plugins = config.get("plugin", [])
    for name in plugin_names if plugin_names else []:
        if name not in PLUGINS:
            raise ValueError(f"Unknown plugin: '{name}'.")

    if command in ["enable", "activate"]:
        [active_plugins.append(name) for name in plugin_names]

    elif command in ["disable", "deactivate"]:
        [active_plugins.remove(name) for name in plugin_names]

    else:
        raise ValueError(
            f"Invalid plugin command: '{command}'. Options: 'list', 'enable' or 'disable'."
        )

    active_plugins = sorted(list(set(active_plugins)))
    config["plugin"] = active_plugins
    config.save()
    print(f"Enabled plugins: {', '.join(active_plugins)}")
