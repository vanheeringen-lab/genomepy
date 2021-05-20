import argparse
from tempfile import TemporaryDirectory

import pytest

import genomepy.argparse_support


def test_argparse_plugin():
    action = genomepy.argparse_support.parse_genome

    with TemporaryDirectory() as tmpdir:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-g", dest="genome", action=action(auto_install=True, genomes_dir=tmpdir)
        )
        args = parser.parse_args(
            ["-g", "ASM2732v1"],
        )
        assert isinstance(args.genome, genomepy.Genome)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parser = argparse.ArgumentParser()
        parser.add_argument("-g", dest="genome", action=action())
        _ = parser.parse_args(["-g", "non_existing"])
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
