import subprocess
import sys

import pytest

from hrdscope.cli import main

COMMANDS = ["labels", "score", "splits", "embed", "stream", "fetch", "work", "tune", "benchmark", "predict"]


@pytest.mark.parametrize("cmd", COMMANDS)
def test_subcommand_help(cmd, capsys):
    with pytest.raises(SystemExit) as exc:
        main([cmd, "--help"])
    assert exc.value.code == 0
    assert cmd in capsys.readouterr().out


def test_score_runs_on_reference_example(tmp_path):
    seg = tmp_path / "x.seg.txt"
    seg.write_text("GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number\n"
                   "a\tchr1\t1\t50000000\t2\t2\t0\na\tchr1\t50000000\t248000000\t2\t1\t1\n")
    out = subprocess.run([sys.executable, "-m", "hrdscope.cli", "score", str(seg)], capture_output=True, text=True)
    assert out.returncode == 0 and "HRD_sum=" in out.stdout


def test_is_tumour_sample():
    from hrdscope.queue import is_tumour_sample

    assert is_tumour_sample("TCGA-04-1331-01A-01-TS1.dcm")
    assert is_tumour_sample("TCGA-13-A5FT-02Z-00-DX1.svs")
    assert not is_tumour_sample("TCGA-61-1734-11A-01-TS1.c55c7e4a.svs")
    assert is_tumour_sample("1009484_171020_ImageActual.svs")
