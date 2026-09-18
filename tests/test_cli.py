import subprocess
import sys

import pytest

from hrdscope.cli import main

COMMANDS = ["labels", "score", "splits", "embed", "stream", "tune", "benchmark", "predict"]


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
