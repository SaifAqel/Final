# tests/test_main_full_trace.py
import sys, runpy
from pathlib import Path

def test_main_full_trace(capsys, monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    main_py = repo_root / "main.py"

    # ensure imports like `from logging_utils import ...` work
    monkeypatch.chdir(repo_root)
    monkeypatch.syspath_prepend(str(repo_root))

    argv = [
        "main.py",
        "--stages", "config/stages.yaml",
        "--streams", "config/streams.yaml",
        "--drum", "config/drum.yaml",
        "--log", "TRACE",
    ]
    monkeypatch.setenv("PYTHONUNBUFFERED", "1")
    monkeypatch.setattr(sys, "argv", argv, raising=False)

    # run the script as __main__
    runpy.run_path(str(main_py), run_name="__main__")

    out, err = capsys.readouterr()
    if out:
        print(out)
    if err:
        print(err)

    assert "Gas T_out:" in out
