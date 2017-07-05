#!/usr/bin/env python

"""
Analysis Tools for CADEE

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import os
import sys
import gzip
import tools
import numpy as np


__author__ = "Beat Amrein, Miha Purg"
__email__ = "beat.amrein@gmail.com, miha.purg@gmail.com"

logger = tools.getLogger()

DEBUG = True
DEBUG = False

NLC = '\n'


try:
    import qscripts
    from qscripts import q_genfeps
    from qscripts import q_mapper
    from qscripts import q_analysemaps
    from qscripts import q_analysedyns
    from qscripts import g_genrelax

except ImportError:
    QSCRIPTS_DIR = os.path.dirname(os.path.realpath(__file__)) + "/../qscripts/"

    if not os.path.isdir(QSCRIPTS_DIR):
        print('FATAL: please copy qscripts to {}'.format(QSCRIPTS_DIR))
        raise Exception('No qscripts: %s', QSCRIPTS_DIR)

    sys.path.insert(1, QSCRIPTS_DIR)

    import q_genfeps       # NOPEP8
    import q_mapper        # NOPEP8
    import q_analysemaps   # NOPEP8
    import q_analysedyns   # NOPEP8
    import q_genrelax      # NOPEP8


class MapSettings(object):
    def __init__(self, mutant, replik, alpha, hij, force_remap=False):
        self.mutant = mutant
        self.replik = replik
        self.alpha = alpha
        self.hij = hij
        if force_remap:
            self.force_remap = True
        else:
            self.force_remap = False


class QScriptsError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class AnalysisInternalError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class AnalysisError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InvalidRankError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def analyse_with_qscripts(logfiles, results):
    """ Analyse Logfiles. Return Results-object.
    :param logfiles: list of logfiles
    :param results: results Object
    :type logfiles: str
    :type results: Results
    :return: results
    :type return: Results
    """

    # TODO: store offdiagonals, store per-file information
    #       (rewrite Results object to make this convenient

    qads = q_analysedyns.QAnalyseDyns(logfiles, timeunit="fs")

    t, pot_e = qads.get_energies("SUM", percent_skip=10).get_columns(["Time", "Potential"])  # NOPEP8
    t, kin_e = qads.get_energies("SUM", percent_skip=10).get_columns(["Time", "Kinetic"])  # NOPEP8
    t, tot_e = qads.get_energies("SUM", percent_skip=10).get_columns(["Time", "Total"])  # NOPEP8

    t, qq_el = qads.get_q_energies("Q-Q", 2, percent_skip=10).get_columns(["Time", "El"])   # NOPEP8
    ts = qads.get_temps(percent_skip=10)
    t_tot = ts.get_columns(["T_tot"])[0]
    t_free = ts.get_columns(["T_free"])[0]
    t_free_solute = ts.get_columns(["T_free_solute"])[0]
    t_free_solvent = ts.get_columns(["T_free_solvent"])[0]

    results.ene(np.mean(kin_e), np.mean(pot_e), np.mean(tot_e))

    results.temp(np.mean(t_tot), np.mean(t_free),
                 np.mean(t_free_solvent), np.mean(t_free_solute))

    return results


def main(mset, eqfil):
    """
    :param mset: mapsettings
    :parame eqfil: map fep preceding eqfil
    """

    def get_number(fname):
        if not isinstance(fname, str):
            return sys.maxint

        fname = os.path.basename(fname)

        num = fname.split("_")[0]

        try:
            return int(num)
        except ValueError:
            return sys.maxint
        except TypeError:
            return sys.maxint

    def reset_mapping():
        """delete .qana.mapped files"""
        logger.info('deleting old .qana.mapped files...')
        for fil in os.listdir('./'):
            if fil.endswith('.qana.mapped'):
                logger.debug('deleting %s', fil)
                os.remove(fil)

    if isinstance(mset, MapSettings):
        if mset.force_remap:
            reset_mapping()
            mset.force_remap = False
    else:
        raise TypeError('Is not a valid MapSettings object', mset)

    max_number = get_number(eqfil)

    if DEBUG:
        logger.debug('max_number %s', max_number)

    qana_fils = []
    for fil in os.listdir('./'):
        if fil.endswith('.qana'):
            if os.path.isfile(fil + '.mapped'):
                if DEBUG:
                    logger.debug('found .mapped %s', fil)
                continue
            qana_fils.append(fil)

    for qana in sorted(qana_fils):

        if eqfil is not None:
            if get_number(qana) > max_number:
                if DEBUG:
                    logger.debug('max_number %s, mynumber %s', max_number,
                                 get_number(qana))
                return False

        fils = open(qana).readline().split()

        try:
            results = map_and_analyse(fils, mset.mutant, mset.replik,
                                      mset.alpha, mset.hij)
            open(qana + ".mapped", 'w').write(str(results.items()))
        except AnalysisInternalError as aie:
            logger.exception('Mapping failed because of reason: %s', aie)
        except AnalysisError as ae:
            logger.exception('Incomlete Energy or Logfiles: %s', ae)
        except QScriptsError as qse:
            logger.exception('QScripts Failed! Reason: %s', qse)
        except q_mapper.QMappingError as qme:
            logger.exception('QMappingError! Reason: %s', qme)
    return True


def map_and_analyse(fils, mutant, replik, alpha, hij):

    cleanup_list = []

    def all_files_exist(files):
        """check if all files in files exist.
        :param files: list of files
        :type files: list
        :return: bool
        """
        if len(files) == 0:
            return False
        for fil in files:
            if not os.path.isfile(fil):
                if os.path.isfile(fil + '.gz'):
                    with gzip.open(fil + '.gz') as gzfil:
                        data = gzfil.read()
                        open(fil, 'w').write(data)
                        cleanup_list.append(fil)
                        continue
                logger.debug('Does not exist %s, %s', fil, fil+'.gz')
                return False
        return True

    def cleanup():
        for tmpfil in cleanup_list:
            os.remove(tmpfil)

    logs = []
    enes = []
    medium = True
    eqfil = None
    for fil in fils:
        logs.append(fil+".log")
        enes.append(fil+".en")

        if fil[3] != "0":
            medium = False

        if "_eq" in fil:
            eqfil = fil

    if eqfil is None:
        logger.warning('Error, no _eq file in file-list')
        raise AnalysisInternalError('Error, no _eq file in file-list')

    # write the energy filenames to a file for QMapper
    open("q_enfiles.list", "w").write("\n".join(enes))

    # parameters for QMapper
    qmap_args = {"hij": hij, "gas_shift": alpha,
                 "bins": 20, "skip": 250,
                 "minpts_per_bin": 20,
                 "temperature": 300,
                 "en_list_fn": "q_enfiles.list",
                 "verbose": False,
                 "mapdirectories": ["./"],
                 "nthreads": 1}

    if medium:
        # medium sized fep
        logger.debug('medium sized fep!')
        qmap_args.update({"minpts_per_bin": 100, "bins": 62})
        fepsize = 'medium'
    elif len(logs) == 11:
        # ultra short fep
        logger.debug('us fep!')
        fepsize = 'us'
    else:
        # tripplet fep
        logger.debug('tripple us fep!')
        qmap_args.update({"minpts_per_bin": 40, "bins": 40})
        fepsize = 'tripple_us'

    if all_files_exist(logs) and all_files_exist(enes):

        # initialize qmapper object
        qmapper = q_mapper.QMapper(**qmap_args)

        # assign a useful name to our results

        results = tools.Results(mutant, replik, eqfil, fepsize)

        results = analyse_with_qscripts(logs, results)

        # map the run and check for failure
        (mapped, failed) = qmapper.q_map()
        if failed:   # list of tuples -> [ (mapdir, error), ... ]
            err_msg = failed[0][1]
            logger.warning(
                'Error while mapping (mutant: %s, '
                'replik: %s, fepsize: %s, eq_file: %s): '
                '%s', mutant, replik, fepsize, eqfil, err_msg)
            cleanup()
            raise QScriptsError('Mapping Failed')

        # analyse the mapped run (qfep output)
        # and check for failure
        qana = q_analysemaps.QAnalyseMaps(mapped)
        if failed:   # list of tuples -> [ (mapdir, error), ... ]
            err_msg = failed[0][1]
            logger.warning(
                'Error while analysing (mutant: %s, '
                'replik: %s, fepsize: %s,eq_file: %s): '
                '%s', mutant, replik, fepsize, eqfil, err_msg)
            cleanup()
            raise QScriptsError('Mapping Analysis Failed')

        try:
            qa1 = qana.get_analysed()[0]
        except IndexError:
            logger.warning(
                'Error while extracting (mutant: %s, '
                'replik: %s, fepsize: %s,eq_file: %s): '
                '%s', mutant, replik, fepsize, eqfil, 'IndexError'
                )
            cleanup()
            raise QScriptsError('Mapping Analysis Failed: IndexError')

        dGa, dG0 = qa1.get_dGa(), qa1.get_dG0()
        # logging
        logger.info('%s %s %s %s dGa: %s dG0: %s',
                    mutant,
                    replik,
                    fepsize,
                    eqfil,
                    dGa, dG0
                    )
        results.dg(dGa, dG0)
        tools.report_results(results)
    else:
        logger.warning(
            'Could not find logfiles!'
            'args %s', logs
            )
        cleanup()
        raise AnalysisError('Could not find all log- and energy files')

    cleanup()
    return results
