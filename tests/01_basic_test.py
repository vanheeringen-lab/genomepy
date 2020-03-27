import genomepy


def test_import():
    # __init__.py
    assert str(genomepy.search).startswith("<function search at")
    assert str(genomepy.Genome) == "<class 'genomepy.genome.Genome'>"
    assert str(genomepy.ProviderBase) == "<class 'genomepy.provider.ProviderBase'>"
    assert genomepy.__author__ == "Simon van Heeringen"


def test_config():
    cfg = genomepy.functions.config
    print(cfg)
    assert len(cfg.keys()) == 3
