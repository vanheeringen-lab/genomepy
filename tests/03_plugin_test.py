import genomepy
import pytest
import norns

config = norns.config("genomepy", default="cfg/default.yaml")


def test_find_plugins():
    genomepy.plugin.find_plugins()
    function_loc = genomepy.plugins.bwa
    assert str(function_loc).endswith("genomepy/genomepy/plugins/bwa.py'>")


def test_convert(name="TestName", expected="test_name"):
    result = genomepy.plugin.convert(name)
    assert result == expected


def test_init_plugins():
    # returns dict of all plugins and their functions
    p = genomepy.plugin.init_plugins()
    expected = ["minimap2", "bowtie2", "gmap", "hisat2", "star", "blacklist", "bwa"]
    assert isinstance(p, dict)
    assert list(p.keys()) == expected


def test_get_active_plugins():
    # returns list of all active plugins
    plugins = genomepy.plugin.get_active_plugins()
    assert isinstance(plugins, list)
    assert config["plugin"] == plugins


def test_activate():
    # test assumes BWA is currently inactive
    plugins = genomepy.plugin.get_active_plugins()
    assert len([plugin for plugin in plugins if "bwa" in str(plugin)]) == 0

    genomepy.plugin.activate("bwa")
    plugins = genomepy.plugin.get_active_plugins()
    assert len([plugin for plugin in plugins if "bwa" in str(plugin)]) == 1


def test_deactivate():
    # test assumes BWA is currently active
    plugins = genomepy.plugin.get_active_plugins()
    assert len([plugin for plugin in plugins if "bwa" in str(plugin)]) == 1

    genomepy.plugin.deactivate("bwa")
    plugins = genomepy.plugin.get_active_plugins()
    assert len([plugin for plugin in plugins if "bwa" in str(plugin)]) == 0


def test_plugin():
    p = genomepy.plugin.Plugin()
    for class_method in [
        "name",
        "activate",
        "deactivate",
        "after_genome_download",
        "get_properties",
    ]:
        assert class_method in dir(p)

    # test name method
    assert p.name() == ""

    # test not implemented method
    with pytest.raises(NotImplementedError):
        p.after_genome_download(genome="", threads=2, force=False)
