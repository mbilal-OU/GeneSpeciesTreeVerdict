from pathlib import Path

import pytest

from genespeciestreeverdict.tutorial import write_tutorial_dataset


@pytest.fixture()
def tutorial_dir(tmp_path: Path) -> Path:
    out = tmp_path / "tutorial"
    write_tutorial_dataset(out)
    return out
