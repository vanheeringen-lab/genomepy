import pytest

import genomepy.plugins


@pytest.fixture
def deactivate_plugins():
    # save originally active plugins
    original_plugins = [p.name for p in genomepy.plugins.get_active_plugins()]
    # deactivate all plugins
    [genomepy.plugins.deactivate(p) for p in original_plugins]

    yield

    # reactivate original plugins
    [genomepy.plugins.activate(p) for p in original_plugins]


def test_plugin():
    p = genomepy.plugins.Plugin()
    assert p.name == "plugin"

    for class_method in [
        "activate",
        "deactivate",
        "after_genome_download",
        "get_properties",
    ]:
        assert class_method in dir(p)

    assert p.active is False
    p.activate()
    assert p.active is True


def test_convert(name="TestName", expected="test_name"):
    result = genomepy.plugins.convert(name)
    assert result == expected


def test_list_plugins():
    plugins = genomepy.plugins.list_plugins()
    expected = ["blacklist", "bowtie2", "bwa", "gmap", "hisat2", "minimap2", "star"]
    assert sorted(plugins) == sorted(expected)


def test_init_plugins():
    p = genomepy.plugins.init_plugins()
    assert isinstance(p, dict)

    plugins = list(p.keys())
    expected = ["blacklist", "bowtie2", "bwa", "gmap", "hisat2", "minimap2", "star"]
    assert sorted(plugins) == sorted(expected)


def test_get_active_plugins(deactivate_plugins):
    _ = deactivate_plugins

    # returns list of all active plugins
    genomepy.plugins.activate("bwa")
    active_plugins = genomepy.plugins.get_active_plugins()
    assert isinstance(active_plugins, list)
    assert ["bwa"] == [p.name for p in active_plugins]


def test_activate(deactivate_plugins):
    _ = deactivate_plugins

    assert len(genomepy.plugins.get_active_plugins()) == 0
    genomepy.plugins.activate("bwa")
    assert len(genomepy.plugins.get_active_plugins()) == 1
    genomepy.plugins.deactivate("bwa")


def test_deactivate(deactivate_plugins):
    _ = deactivate_plugins

    genomepy.plugins.activate("bwa")
    assert len(genomepy.plugins.get_active_plugins()) == 1
    genomepy.plugins.deactivate("bwa")
    assert len(genomepy.plugins.get_active_plugins()) == 0


def test_manage_plugins(capsys):
    genomepy.plugins.manage_plugins("enable", ["blacklist"])
    genomepy.plugins.manage_plugins("list")
    captured = capsys.readouterr().out.strip().split("\n")
    assert captured[2].startswith("blacklist")
    assert captured[2].endswith("*")

    genomepy.plugins.manage_plugins("disable", ["blacklist"])
    genomepy.plugins.manage_plugins("list")
    captured = capsys.readouterr().out.strip().split("\n")
    assert captured[2].startswith("blacklist")
    assert not captured[2].endswith("*")

    with pytest.raises(ValueError):
        genomepy.plugins.manage_plugins("blurp")
