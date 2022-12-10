#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `structure_step` package."""

import pytest  # noqa: F401
import structure_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = structure_step.Structure()
    assert str(type(result)) == "<class 'structure_step.structure.Structure'>"
