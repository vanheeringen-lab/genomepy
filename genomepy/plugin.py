import os
import sys
import re

import norns

config = norns.config("genomepy", default="cfg/default.yaml")

class Plugin(object):
    """Plugin base class.
    """
    active = False

    def name(self):
        n = type(self).__name__.replace("Plugin", "")
        return convert(n)

    def activate(self):
        self.active = True

    def deactivate(self):
        self.active = False

    def after_genome_download(self, genome):
        raise NotImplementedError("plugin should implement this method")

    def get_properties(self, genome):
        raise NotImplementedError("plugin should implement this method")

def find_plugins():
    """Locate and initialize all available plugins.
    """ 
    plugin_dir = os.path.dirname(os.path.realpath(__file__))
    plugin_dir = os.path.join(plugin_dir, "plugins")
    plugin_files = [x[:-3] for x in os.listdir(plugin_dir) if x.endswith(".py")]
    for plugin in plugin_files:
        __import__("genomepy.plugins.{}".format(plugin))

def convert(name):
    """Convert CamelCase to underscore

    Parameters
    ----------
    name : str
        Camelcase string

    Returns
    -------
    name : str
        Converted name
    """ 
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

def init_plugins():
    """Return dictionary of available plugins

    Returns
    -------
    plugins : dictionary
        key is plugin name, value Plugin object
    """ 
    find_plugins()
    d = {}
    for c in Plugin.__subclasses__():
        ins = c()
    
        if ins.name() in config.get("plugin", []):
            ins.activate()
        
        d[ins.name()] = ins
    
    return d

def activate(name):
    """Activate plugin.

    Parameters
    ----------
    name : str
        Plugin name.
    """
    if name in plugins:
        plugins[name].activate()
    else:
        raise Exception("plugin {} not found".format(name))

def deactivate(name):
    """Deactivate plugin.

    Parameters
    ----------
    name : str
        Plugin name.
    """
    if name in plugins:
        plugins[name].deactivate()
    else:
        raise Exception("plugin {} not found".format(name))

def get_active_plugins():
    """Returns all active plugin instances.
    """ 
    return [inst for name, inst in plugins.items() if inst.active]

plugins = init_plugins()
