# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
Tests for tsId disambiguation when importing tilt-series from mdoc files
with duplicate basenames across different subdirectories.
"""
import os
import unittest

from tomo.protocols.protocol_ts_import import (
    _disambiguateTsIds,
    _sanitizeTsIdComponent,
)


class TestSanitizeTsIdComponent(unittest.TestCase):
    """Unit tests for _sanitizeTsIdComponent."""

    def test_no_special_chars(self):
        self.assertEqual(_sanitizeTsIdComponent('dir1'), 'dir1')

    def test_hyphen_replaced(self):
        self.assertEqual(_sanitizeTsIdComponent('my-dir'), 'my_dir')

    def test_dot_removed(self):
        self.assertEqual(_sanitizeTsIdComponent('my.dir'), 'mydir')

    def test_brackets_removed(self):
        self.assertEqual(_sanitizeTsIdComponent('[dir]'), 'dir')

    def test_double_underscore_collapsed(self):
        self.assertEqual(_sanitizeTsIdComponent('a__b'), 'a_b')

    def test_combined_special_chars(self):
        self.assertEqual(_sanitizeTsIdComponent('my-dir.2[a]'), 'my_dir2a')


class TestDisambiguateTsIds(unittest.TestCase):
    """Unit tests for _disambiguateTsIds."""

    def test_no_duplicates_unchanged(self):
        """When all tsIds are unique, they should be returned as-is."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir1/position1.mdoc'),
            ('position2', '/data/acquisition/dir1/position2.mdoc'),
            ('position3', '/data/acquisition/dir2/position3.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(result, ['position1', 'position2', 'position3'])

    def test_duplicate_basenames_different_dirs(self):
        """Core scenario: same basename in different subdirectories produces
        unique tsIds by prepending the sanitized relative directory."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir1/position1.mdoc'),
            ('position2', '/data/acquisition/dir1/position2.mdoc'),
            ('position1', '/data/acquisition/dir2/position1.mdoc'),
            ('position2', '/data/acquisition/dir2/position2.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 4)
        # All ids must be unique
        self.assertEqual(len(set(result)), 4)
        # Disambiguated ids should include directory prefix
        self.assertIn('dir1_position1', result)
        self.assertIn('dir2_position1', result)
        self.assertIn('dir1_position2', result)
        self.assertIn('dir2_position2', result)

    def test_mixed_unique_and_duplicate(self):
        """Only duplicated tsIds get disambiguated; unique ones stay as-is."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir1/position1.mdoc'),
            ('position1', '/data/acquisition/dir2/position1.mdoc'),
            ('unique_ts', '/data/acquisition/dir3/unique_ts.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 3)
        self.assertEqual(len(set(result)), 3)
        self.assertIn('dir1_position1', result)
        self.assertIn('dir2_position1', result)
        self.assertEqual(result[2], 'unique_ts')

    def test_nested_subdirectories(self):
        """Nested directories are flattened with underscores."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/session1/grid1/position1.mdoc'),
            ('position1', '/data/acquisition/session1/grid2/position1.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(set(result)), 2)
        self.assertIn('session1_grid1_position1', result)
        self.assertIn('session1_grid2_position1', result)

    def test_directory_with_special_chars(self):
        """Special characters in directory names are sanitized."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir-1/position1.mdoc'),
            ('position1', '/data/acquisition/dir-2/position1.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(set(result)), 2)
        self.assertIn('dir_1_position1', result)
        self.assertIn('dir_2_position1', result)

    def test_digit_starting_directory(self):
        """If the disambiguated tsId starts with a digit, TS_ prefix is added."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/1dir/position1.mdoc'),
            ('position1', '/data/acquisition/2dir/position1.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(set(result)), 2)
        # Both should get TS_ prefix since they start with a digit
        self.assertIn('TS_1dir_position1', result)
        self.assertIn('TS_2dir_position1', result)

    def test_empty_list(self):
        """Empty input returns empty output."""
        result = _disambiguateTsIds([], '/data')
        self.assertEqual(result, [])

    def test_single_entry(self):
        """A single entry cannot be a duplicate."""
        fpath = '/data/acquisition'
        pairs = [('position1', '/data/acquisition/dir1/position1.mdoc')]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(result, ['position1'])

    def test_three_duplicates(self):
        """Three tilt-series with the same basename in different directories."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir1/position1.mdoc'),
            ('position1', '/data/acquisition/dir2/position1.mdoc'),
            ('position1', '/data/acquisition/dir3/position1.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 3)
        self.assertEqual(len(set(result)), 3)
        self.assertIn('dir1_position1', result)
        self.assertIn('dir2_position1', result)
        self.assertIn('dir3_position1', result)

    def test_deterministic_ordering(self):
        """Results preserve the input order."""
        fpath = '/data/acquisition'
        pairs = [
            ('position1', '/data/acquisition/dir2/position1.mdoc'),
            ('position1', '/data/acquisition/dir1/position1.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(result[0], 'dir2_position1')
        self.assertEqual(result[1], 'dir1_position1')

    def test_ts_prefix_preserved_in_original(self):
        """If the original tsId already has TS_ prefix (from normalization),
        disambiguation still works correctly."""
        fpath = '/data/acquisition'
        pairs = [
            ('TS_001', '/data/acquisition/dir1/001.mdoc'),
            ('TS_001', '/data/acquisition/dir2/001.mdoc'),
        ]
        result = _disambiguateTsIds(pairs, fpath)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(set(result)), 2)
        self.assertIn('dir1_TS_001', result)
        self.assertIn('dir2_TS_001', result)


if __name__ == '__main__':
    unittest.main()
