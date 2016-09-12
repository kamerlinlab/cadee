#!/usr/bin/env python
"""
This are unittests for tools.py

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import unittest
import tools as tools
import os

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

class MyToolsTests(unittest.TestCase):
    waterpdb = """REMARK THIS IS A TEST WITH WATER ONLY
ATOM   5341  O   HOH   466      30.350 117.099  41.362
ATOM   5342  H1  HOH   466      30.351 116.208  41.710
ATOM   5343  H2  HOH   466      30.035 117.006  40.463
ATOM   5344  O   HOH   467      30.350 120.201  22.750
ATOM   5345  H1  HOH   467      30.744 119.572  23.353
ATOM   5346  H2  HOH   467      29.482 120.374  23.115
ATOM   5347  O   H2O   468      30.350 120.201  32.056
ATOM   5348  H1  H2O   468      30.804 119.606  31.460
ATOM   5349  H2  H2O   468      29.961 119.629  32.717
ATOM   5350  O   WAT   469      30.350 120.201  35.158
ATOM   5351  H1  WAT   469      30.440 120.765  35.926
ATOM   5352  H2  WAT   469      29.625 120.584  34.664
ATOM   5353  O   HOH   470      33.452  89.181  35.158
ATOM   5354  H1  HOH   470      33.864  89.527  35.950
"""

    def test_isint(self):
        self.assertFalse(tools.isint('a'))
        self.assertFalse(tools.isint('0.3'))
        self.assertTrue(tools.isint('3'))
        self.assertTrue(tools.isint('-1'))

    def test_isnum(self):
        self.assertTrue(tools.isnum('0.1'))
        self.assertTrue(tools.isnum('-0.1'))
        self.assertTrue(tools.isnum('0.0e10'))
        self.assertTrue(tools.isnum('5'))
        self.assertFalse(tools.isnum('inf'))
        self.assertFalse(tools.isnum('nan'))
        self.assertFalse(tools.isnum('ab'))

    def test_log_exists(self):
        import os
        fil = '.test.logfile.temporary'
        with open(fil, 'w') as f:
            f.write('\n')
        self.assertTrue(tools.bool_log_exists(fil))
        os.remove(fil)
        self.assertFalse(tools.bool_log_exists(fil))

    def test_pdb_water_only(self):
        import os
        fil = '.test.logfile.temporary.pdb'
        with open(fil, 'w') as f:
            f.write(MyToolsTests.waterpdb)
        self.assertTrue(tools.is_pdb_water_only(fil))
        with open(fil, 'a') as f:
            f.write("""ATOM   5342  O   ASP   467      30.350 117.099  41.362
                    """)
        self.assertFalse(tools.is_pdb_water_only(fil))
        os.remove(fil)

    def test_euklid_dist(self):
        dot0 = [0, 0, 0]
        dot1 = [0, 0, 1]
        dot2 = [0, 0, 2]
        self.assertEqual(tools.euklid_dist(dot0, dot1), 1)
        self.assertEqual(tools.euklid_dist(dot0, dot2), 2)
        self.assertEqual(tools.euklid_dist(dot1, dot2), 1)

    def test_get_atomnumber(self):
        fil = '.test.logfile.temporary.pdb'
        with open(fil, 'w') as f:
            f.write(MyToolsTests.waterpdb)
        self.assertEqual(tools.get_atomnumber(fil, 'O', 'HOH', 470), 5353)
        with open(fil, 'a') as f:
            f.write("""ATOM   5362  O   ASP   467      30.350 117.099  41.362
""")
            f.write("""HETATM 5363  O   ASP   468      30.350 117.099  41.362
""")
            f.write("""ATOM   5364  O   ASP   469      30.350 117.099  41.362
""")
        self.assertEqual(tools.get_atomnumber(fil, 'O', 'ASP', 467), 5362)
        self.assertEqual(tools.get_atomnumber(fil, 'O', 'ASP', 468), 5363)

        # INEXISTENT RECORD
        with self.assertRaisesRegexp(Exception, "Atom.*not found in.*"):
            tools.get_atomnumber(fil, 'O', 'HIS', 466)

        # NOT ATOM RECORD
        with open(fil, 'a') as f:
            f.write("""TER    5362  O   ASP   666
""")
        with self.assertRaisesRegexp(Exception, "Atom.*not found in.*"):
            tools.get_atomnumber(fil, 'O', 'ASP', 666)

        os.remove(fil)

    def test_rename_pdb_res(self):
        residue = ['ATOM   5341  O   HOH   466      30.350 117.099  41.362',
                   'ATOM   5342  H1  HOH   466      30.351 116.208  41.710',
                   'ATOM   5343  H2  HOH   466      30.035 117.006  40.463']
        renres = ['ATOM   5341  O   H2O   466      30.350 117.099  41.362',
                  'ATOM   5342  H1  H2O   466      30.351 116.208  41.710',
                  'ATOM   5343  H2  H2O   466      30.035 117.006  40.463']

        self.assertEqual(tools.rename_pdb_res(residue, 'H2O'), renres)

        with self.assertRaisesRegexp(Exception, 'PDB-Residue must be exactly 3 characters.'):  # NOPEP8
            tools.rename_pdb_res(residue, 'WATER')

    def test_check_qprep_pdb(self):
        line1 = 'ATOM   5342  H1  HOH   466      30.351 116.208  41.710'
        line2 = 'HETATM 5343  H1  HOH   466      30.351 116.208  41.710'
        line3 = 'TER    5344  H1  HOH   466      30.351 116.208  41.710'
        line4 = 'GAP    5344  H1  HOH   466      30.351 116.208  41.710'
        self.assertTrue(tools.check_qprep_pdb(line1))
        self.assertTrue(tools.check_qprep_pdb(line2))
        self.assertFalse(tools.check_qprep_pdb(line3))
        self.assertFalse(tools.check_qprep_pdb(line4))

        line5 = 'ATOM   5342  H1  HOH   466      30.351 116.208  41.710  '
        with self.assertRaisesRegexp(Exception, "Not a Qprep5-ed pdbfile.*"):
            tools.check_qprep_pdb(line5)
        line6 = 'HETATM 5343  H1  HOH   466      30.351 116.208  41.710 blabalba'  # NOPEP8
        with self.assertRaisesRegexp(Exception, "Not a Qprep5-ed pdbfile.*"):
            tools.check_qprep_pdb(line6)

unittest.main()
