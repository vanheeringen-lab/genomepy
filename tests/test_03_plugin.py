import pytest

import genomepy.plugin


@pytest.fixture
def deactivate_plugins():
    # save originally active plugins
    original_plugins = [p.name() for p in genomepy.plugin.get_active_plugins()]
    # deactivate all plugins
    [genomepy.plugin.deactivate(p) for p in original_plugins]

    yield

    # reactivate original plugins
    [genomepy.plugin.activate(p) for p in original_plugins]


def test_find_plugins():
    genomepy.plugin.find_plugins()
    plugin_dict = genomepy.plugin.plugins
    assert "genomepy.plugins.bwa.BwaPlugin object at" in str(plugin_dict["bwa"])


def test_convert(name="TestName", expected="test_name"):
    result = genomepy.plugin.convert(name)
    assert result == expected


def test_init_plugins():
    # returns dict of all plugins and their functions
    p = genomepy.plugin.init_plugins()
    result = sorted(list(p.keys()))
    expected = ["blacklist", "bowtie2", "bwa", "gmap", "hisat2", "minimap2", "star"]
    assert isinstance(p, dict)
    assert result == expected


def test_get_active_plugins(deactivate_plugins):
    _ = deactivate_plugins

    # returns list of all active plugins
    genomepy.plugin.activate("bwa")
    active_plugins = genomepy.plugin.get_active_plugins()
    assert isinstance(active_plugins, list)
    assert ["bwa"] == [p.name() for p in active_plugins]


def test_activate(deactivate_plugins):
    _ = deactivate_plugins

    assert len(genomepy.plugin.get_active_plugins()) == 0
    genomepy.plugin.activate("bwa")
    assert len(genomepy.plugin.get_active_plugins()) == 1


def test_deactivate(deactivate_plugins):
    _ = deactivate_plugins

    genomepy.plugin.activate("bwa")
    assert len(genomepy.plugin.get_active_plugins()) == 1
    genomepy.plugin.deactivate("bwa")
    assert len(genomepy.plugin.get_active_plugins()) == 0


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
