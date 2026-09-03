"""Tests for retired OxygenPriorityAssignment stub."""

import pytest

from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment


def test_oxygen_priority_assignment_raises():
    with pytest.raises(NotImplementedError, match="no longer supported"):
        OxygenPriorityAssignment(None)
