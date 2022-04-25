#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass
from shutil import copy
from filecmp import cmp

import pytest

from itaxotools.sequence_bouncer import SequenceBouncer

TEST_DATA_DIR = Path(__file__).parent


@dataclass
class AcceptanceTest:
    name: str
    input: str


test_data = [
    AcceptanceTest('simple', 'simple.fasta'),
    ]


@pytest.mark.parametrize("test", test_data)
def test_acceptance(test: AcceptanceTest, tmp_path: Path) -> None:

    test_path = TEST_DATA_DIR / test.name
    test_input = test_path / test.input
    tmp_input = tmp_path / test.input
    copy(str(test_input), str(tmp_input))

    SequenceBouncer(str(tmp_input))()

    for test_file in test_path.iterdir():
        tmp_file = tmp_path / test_file.name
        assert tmp_file.exists
        if test_file.suffix in ['.csv', '.fasta']:
            assert cmp(test_file, tmp_file)
