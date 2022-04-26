#!/usr/bin/env python3

from .SequenceBouncer import (
    SequenceBouncer,
    InputSequence,
    InvalidInput,
    AllColumnsRemovedAsGaps,
    version
    )

__version__ = f'{version}.1'
__all__ = [
    'SequenceBouncer',
    'InputSequence',
    'InvalidInput',
    'AllColumnsRemovedAsGaps'
    ]
