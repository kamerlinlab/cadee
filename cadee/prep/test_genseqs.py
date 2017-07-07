#!/usr/bin/env python

"""
This are the unittests for genseq

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import unittest
import genseqs as genseqs

import logging

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.genseqs')


class MyGenseqsTests(unittest.TestCase):
    ALL = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
           'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def test_genseq2_empty(self):
        first = [MyGenseqsTests.ALL]
        second = genseqs.genseq2(MyGenseqsTests.ALL, [])
        self.assertEqual(first, second)

    def test_genseq2_mutate0(self):
        """try: mutate resid0 (DOES NOT EXISTS)"""
        with self.assertRaises(ValueError):
            genseqs.genseq2(MyGenseqsTests.ALL, [(0, list('ARN'))])

    def test_genseq2_mutate1(self):
        """mutate resid1 to R"""
        first = MyGenseqsTests.ALL[:]
        first[0] = 'R'
        first = [MyGenseqsTests.ALL, first]
        second = genseqs.genseq2(MyGenseqsTests.ALL, [(1, ['R'])])
        self.assertListEqual(first, second)

    def test_genseq2_mutate2(self):
        """only single-letter amino acids are accepted """
        with self.assertRaises(ValueError):
            genseqs.genseq2(MyGenseqsTests.ALL, [(1, ['RS'])])

    def test_genseq2_mutate3(self):
        """introduce 2 mutations on 1 resid """
        a = MyGenseqsTests.ALL[:]
        b = a[:]
        a[0] = 'R'
        b[0] = 'S'
        first = [MyGenseqsTests.ALL, a, b]
        second = genseqs.genseq2(MyGenseqsTests.ALL, [(1, ['R', 'S'])])
        self.assertListEqual(first, second)

    def test_genseq2_mutate4(self):
        """mutate 2 residues on 1 positions"""
        wt = MyGenseqsTests.ALL[0:3]
        a = wt[:]
        b = a[:]
        a[0] = 'T'
        b[1] = 'S'
        c = a[:]
        c[1] = 'S'
        first = [wt, b, a, c]
        second = genseqs.genseq2(wt, [(1, ['T']), (2, ['S'])])
        self.assertListEqual(first, second)

    def test_genseq2_mutate5(self):
        """mutate 2 residues on 2 positions"""
        self.maxDiff = None
        wt = list('ARN')
        a1 = list('ASN')
        a2 = list('AVN')
        a3 = list('TRN')
        a4 = list('QRN')
        a5 = list('TSN')
        a6 = list('QSN')
        a7 = list('TVN')
        a8 = list('QVN')

        first = [wt, a1, a2, a3, a4, a5, a6, a7, a8]
        second = genseqs.genseq2(wt, [(1, ['T', 'Q']), (2, ['S', 'V'])])

        self.assertListEqual(first, second)

    def test_genseq2_mutate6(self):
        """mutate 2 residues on 2 positions, and 1 residue on 1 position"""
        self.maxDiff = None
        wt = list('ARN')
        a1 = list('ASN')
        a2 = list('AVN')
        a3 = list('TRN')
        a4 = list('QRN')
        a5 = list('TSN')
        a6 = list('QSN')
        a7 = list('TVN')
        a8 = list('QVN')

        bt = list('ARW')
        b1 = list('ASW')
        b2 = list('AVW')
        b3 = list('TRW')
        b4 = list('QRW')
        b5 = list('TSW')
        b6 = list('QSW')
        b7 = list('TVW')
        b8 = list('QVW')

        first = [wt, a1, a2, a3, a4, a5, a6, a7, a8,
                 bt, b1, b2, b3, b4, b5, b6, b7, b8]

        second = genseqs.genseq2(wt, [(1, ['T', 'Q']), (2, ['S', 'V']), (3, ['W'])])

        self.assertListEqual(sorted(first), sorted(second))

    def test_genseq2_mutate7(self):
        """mutate 3 residues on 2 positions"""
        self.maxDiff = None
        wt = list('ARN')
        a1 = list('ASN')
        a2 = list('AVN')
        a3 = list('TRN')
        a4 = list('QRN')
        a5 = list('TSN')
        a6 = list('QSN')
        a7 = list('TVN')
        a8 = list('QVN')

        bt = list('ARW')
        b1 = list('ASW')
        b2 = list('AVW')
        b3 = list('TRW')
        b4 = list('QRW')
        b5 = list('TSW')
        b6 = list('QSW')
        b7 = list('TVW')
        b8 = list('QVW')

        ct = list('ARL')
        c1 = list('ASL')
        c2 = list('AVL')
        c3 = list('TRL')
        c4 = list('QRL')
        c5 = list('TSL')
        c6 = list('QSL')
        c7 = list('TVL')
        c8 = list('QVL')

        first = [wt, a1, a2, a3, a4, a5, a6, a7, a8,
                 bt, b1, b2, b3, b4, b5, b6, b7, b8,
                 ct, c1, c2, c3, c4, c5, c6, c7, c8]

        second = genseqs.genseq2(wt, [(1, ['T', 'Q']), (2, ['S', 'V']), (3, ['W', 'L'])])

        self.assertListEqual(sorted(first), sorted(second))

    def test_genseq2_mutate8(self):
        """mutate 3 residues to 3 positions"""
        wt = list('ARN')
        nummut = genseqs.genseq2(wt, [(1, ['T', 'Q', 'P']),
                                      (2, ['S', 'V', 'E']),
                                      (3, ['W', 'L', 'D'])])
        self.assertEqual(len(nummut), 64)

    def test_genseq2_mutate9(self):
        """mutate 3 residues to 3 positions"""
        wt = MyGenseqsTests.ALL
        nummut = genseqs.genseq2(wt, [(1, ['T', 'Q', 'P']),
                                      (9, ['S', 'V', 'E']),
                                      (20, ['W', 'L', 'D'])])
        self.assertEqual(len(nummut), 64)

    def test_genseq2_mutate10(self):
        """mutate 3 residues to 20 positions"""
        wt = MyGenseqsTests.ALL
        all = MyGenseqsTests.ALL
        nummut = genseqs.genseq2(wt, [(1, all),
                                      (9, all),
                                      (20, all)])
        self.assertEqual(len(nummut), 8000)

    def test_genseq2_mutate11(self):
        """mutate 4 residues to 10 positions"""
        wt = MyGenseqsTests.ALL
        ten = MyGenseqsTests.ALL[:10]
        nummut = genseqs.genseq2(wt, [(1, ten),
                                      (3, ten),
                                      (5, ten),
                                      (9, ten)], keepdupes=True)
        self.assertEqual(len(nummut), 11*11*11*11)

unittest.main()
