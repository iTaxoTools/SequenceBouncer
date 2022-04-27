#!/usr/bin/env python3

from pathlib import Path
from typing import Dict
from dataclasses import dataclass
from shutil import copy
from filecmp import cmp

from Bio import SeqIO

import pytest

from itaxotools.sequence_bouncer import SequenceBouncer, InputSequence

TEST_DATA_DIR = Path(__file__).parent


@dataclass
class AcceptanceTest:
    name: str
    input: str
    format: str
    accepted: Dict[str, bool]


test_data = [
    AcceptanceTest('simple', 'simple.fasta', 'fasta', {
        'Sample1': True,
        'Sample2': True,
        'Sample3': True,
        'SampleX': False,
    }),
    ]


@pytest.mark.parametrize("test", test_data)
def test_acceptance(test: AcceptanceTest, tmp_path: Path) -> None:

    test_path = TEST_DATA_DIR / test.name
    test_input = test_path / test.input
    tmp_input = tmp_path / test.input
    copy(str(test_input), str(tmp_input))

    accepted = SequenceBouncer(str(tmp_input))()

    for test_file in test_path.iterdir():
        tmp_file = tmp_path / test_file.name
        assert tmp_file.exists()
        if test_file.suffix in ['.csv', '.fasta']:
            assert cmp(test_file, tmp_file)

    for sample, verdict in accepted.items():
        assert test.accepted[sample] == verdict


@pytest.mark.parametrize("test", test_data)
def test_acceptance_none_path(test: AcceptanceTest, tmp_path: Path) -> None:

    test_path = TEST_DATA_DIR / test.name
    test_input = test_path / test.input
    tmp_input = tmp_path / test.input
    copy(str(test_input), str(tmp_input))

    accepted = SequenceBouncer(Path(tmp_input), write_none=True)()

    for tmp_file in tmp_path.iterdir():
        assert tmp_file.name == test.input

    for sample, verdict in accepted.items():
        assert test.accepted[sample] == verdict


@pytest.mark.parametrize("test", test_data)
def test_acceptance_none_file(test: AcceptanceTest, tmp_path: Path) -> None:

    test_path = TEST_DATA_DIR / test.name
    test_input = test_path / test.input
    tmp_input = tmp_path / test.input
    copy(str(test_input), str(tmp_input))

    with open(tmp_input) as file:
        accepted = SequenceBouncer(file, write_none=True)()

    for tmp_file in tmp_path.iterdir():
        assert tmp_file.name == test.input

    for sample, verdict in accepted.items():
        assert test.accepted[sample] == verdict


@pytest.mark.parametrize("test", test_data)
def test_acceptance_none_input(test: AcceptanceTest, tmp_path: Path) -> None:

    test_path = TEST_DATA_DIR / test.name
    test_input = test_path / test.input
    tmp_input = tmp_path / test.input
    copy(str(test_input), str(tmp_input))

    iterator = SeqIO.parse(tmp_input, test.format)
    input = InputSequence('gene', iterator)
    accepted = SequenceBouncer(input, write_none=True)()

    for tmp_file in tmp_path.iterdir():
        assert tmp_file.name == test.input

    for sample, verdict in accepted.items():
        assert test.accepted[sample] == verdict
