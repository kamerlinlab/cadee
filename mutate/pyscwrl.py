#!/usr/bin/env python

"""Run SCWRL4 for mutagenesis

This module is interfacing with the Scwrl4 executable.

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import os
import shutil
import time

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.pyscwrl')

NLC = '\n'


def run_scwrl(out_pdb, out_log, in_pdb, in_seq):
    """Run SCWRL4. Use 'Scwrl4' - executable from $PATH"""
    import executables as exe

    scwrl_exe = exe.which('Scwrl4')
    if scwrl_exe is None:
        raise Exception("FATAL: 'Scwrl4' - executable not found.")

    scwrl_line = '{SCWRL} -i {SCWRL_PDB} -o {OUT} -s {SEQ} > {LOG}'
    cmd = scwrl_line.format(SCWRL=scwrl_exe, SCWRL_PDB=in_pdb,
                            OUT=out_pdb, LOG=out_log, SEQ=in_seq)
    os.system(cmd)

    if not os.path.isfile(out_pdb):
        raise Exception("FATAL: 'Scwrl4' - executable did not create output!", cmd)

    # TODO: errorhandling


def post_process_scwrlpdb(resids, wtpdb, scwrlpdb, newpdb):
    """Extract resids from scwrlpdb and replace in wtpdb. Save as newpdb."""
    # delete some contents of newpdb
    open(newpdb, 'w').write('')
    import tools as tools

    def get_residue_from_scwrlpdb(resid, scwrlpdb):
        """
        23 - 26        Integer         Residue sequence number."""
        residue = []
        for line in tools.read_pdbatoms(scwrlpdb):
            if int(line[22:26]) == int(resid):
                residue.append(line.rstrip())
        return residue

    def write_out(line):
        """Append line to outputfile"""
        if 'TER' == line[:3] or not tools.is_hydrogen(line):
            open(newpdb, 'a').write(line.rstrip() + NLC)

    def find_his_protonation(scwrlresidue):
        """SPECIAL CASE FOR HISTIDINE"""
        protons = []
        his = None    # make debugging easy, if HIS not set later
        for scwrlline in scwrlresidue:
            if scwrlline[12:16].strip() == 'HE2':
                protons.append('HE2')
            if scwrlline[12:16].strip() == 'HD1':
                protons.append('HD1')
        if len(protons) == 2:
            if ((protons[0] == 'HE2' or protons[1] == 'HE2') and
                    (protons[0] == 'HD1' or protons[1] == 'HD1')):
                his = 'HIP'
            else:
                raise 'WTF 3'
        elif len(protons) == 1:
            if protons[0] == 'HE2':
                his = 'HIE'
            elif protons[0] == 'HD1':
                his = 'HID'
            else:
                raise 'WTF 1'
        else:
            raise 'WTF 2'

        if his is None:
            raise 'WTF 3'

        return his

    for line in open(wtpdb):
        skip = False
        # slow way
        for resid in resids:
            if int(line[22:26]) == int(resid):
                skip = True
                # pdb info:
                # 13 - 16        Atom            Atom name.
                # 18 - 20        Residue name    Residue name.
                if line[12:16].strip() == 'CA':
                    scwrlresidue = get_residue_from_scwrlpdb(resid, scwrlpdb)
                    if scwrlresidue[0][17:20].strip() == 'HIS':
                        his = find_his_protonation(scwrlresidue)
                        scwrlresidue = tools.rename_pdb_res(scwrlresidue, his)
                    for scwrlline in scwrlresidue:
                        write_out(scwrlline)
        if skip:
            continue
        write_out(line)


def babel_pdb_for_scwrl(wtpdb, proper_pdb='proper.pdb',
                        proper_fasta='proper.fasta', temp_pdb='temp.pdb',
                        sprss_babel=True):
    """ convert q-pdb to scwrl-compatible pdb """

    if sprss_babel:
        sprssbabel = '2>&1 | egrep -v "1 molecule converted|audit log messages"'  # NOPEP8
    else:
        sprssbabel = ''

    while os.path.exists(temp_pdb):
        temp_pdb += '.pdb'

    temp_pdb = '"' + temp_pdb + '"'

    import config
    import executables

    # locate babel-exe
    babel_exe = executables.exe.which('babel')
    if babel_exe is None:
        raise Exception("FATAL: babel executable not found")

    # prepare list of renamed residues
    resren = ''
    for rren in config.RESRENAME:
        resren += 'sed -i "s/{0}/{1}/g" $f'.format(rren[0], rren[1]) + NLC
        if not rren[1] in config.NATURAL_AA:
            raise Exception('Can not replace '+rren[0]+' with' + rren[1]+' which is not recognized by scwrl4!)')  # NOPEP8

    if len(config.NATURAL_AA) < 1:
        raise Exception('Config is not OK. Check config.RESSCWRL4')

    allnat = 'egrep "^(ATOM|HETATM).*('

    for amino in config.NATURAL_AA:
        allnat += amino + '|'

    # kill the last |
    if allnat[-1] == '|':
        allnat = allnat[:-1]
    else:
        raise Exception('WTF')

    allnat += ')" ' + proper_pdb + ' > $f '

    cmd = """

f={TEMP_PDB}

{BABEL} -d ---errorlevel 0 {WT} $f  {SPRSSBABEL}

# convert all protonated residues into standard residues
{RESREN}

# kill protons
sed -i '/H  $/d' $f

# remove connect-records
sed -i '/CONECT/d' $f

mv $f {PROPER_PDB}

# kill everything that is not an natural amino acid
{ALLNAT}

{BABEL} ---errorlevel 1 $f {PROPER_PDB}  {SPRSSBABEL}

/bin/rm $f

{BABEL} ---errorlevel 1 {PROPER_PDB} {PROPER_FASTA} {SPRSSBABEL}

exit 0

""".format(WT=wtpdb, PROPER_PDB=proper_pdb, PROPER_FASTA=proper_fasta,
           TEMP_PDB=temp_pdb, BABEL=babel_exe, SPRSSBABEL=sprssbabel,
           RESREN=resren, ALLNAT=allnat)
    logger.debug(cmd)
    os.system(cmd)


def kill_residue(pdbfile, resnum):
    """Delete Solvent-Residue with resnum resnum in pdbfile
    EDITS pdbfile on disk!
    """

    temp = pdbfile + ".temp"
    ofil = open(temp, 'w')

    import config as config

    for line in open(pdbfile):
        try:
            if int(line[22:26]) == int(resnum):
                resname = line[17:20].strip()
                if resname in config.SOLVENT_MOLECULES:
                    logger.debug('remove: line', line)
                    continue
                else:
                    logger.debug('can not remove %s is not a solvent molecule',
                                 resname)
        except UserWarning:
            pass
        ofil.write(line)
    shutil.move(temp, pdbfile)


def scwrl(wtpdb, wtfep, qprep5inp, mutcode, immutable_resids=[]):
    """
    """
    mutants=[]
    parts = mutcode.split(':')
    resid = int(parts[0])
    if resid in immutable_resids:
        logger.fatal("Fatal: Can not mutate %s, atoms of residue are in immutable_residues!", resid)
        return
    lname = parts[1].upper()
    from config import SatLibs
    aalib = SatLibs.get_lib(lname)
    logger.debug('Will mutate residue %s to %s', resid, aalib)

    mutants.append((resid, aalib))

    from genseqs import get_fasta
    fasta = get_fasta(wtpdb)

    logger.info('wt  sequence:', fasta)
    logger.info('new sequence:', mutants)

    scwrl2(wtpdb, fasta, wtfep, qprep5inp, immutable_resids)


def scwrl2(wtpdb, seq, wtfep, qprep5inp, immutable_resids=[]):
    """Use scrwl to mutate wtpdb into pdb with sequence seq

    :param wtpdb: wild-type pdb used as reference
    :param wtfep: wild-type fep used as reference
    :param qprep5inp: qprep5-inputfile used to generate topology
    :param immutable_resids: residues that will not be touched.
    :type wtpdb: str
    :type wtfep: str
    :type qprep5inp: str
    :type immutable_resids: list[int]

    TODO: -replace babel call with python code

    Assumptions: babel is in $PATH
                 scwrl is can be located by executables.exe

    Output:      write mutant.pdb, mutant.fep, mutant.top in cwd
    """

    import tools
    immutable_resids.extend(tools.get_fep_resids(wtpdb, wtfep))

    logger.debug('wtpdb&fep: %s %s, qprepinp %s, seq: %s', wtpdb, wtfep, qprep5inp, seq)

    def get_resids(seq):
        resids = []
        for i, char in enumerate(seq):
            if char.istitle():
                resids.append(i+1)
        return resids

    start = time.time()

    resids = get_resids(seq)
    log_s4_out = "scwrl.log"
    pdb_s4_out = "scwrl.pdb"
    babeled_pdb = "babeled.pdb"
    seq_s4_in = 'scwrl.seq'
    new4q = "new4q.pdb"

    logger.debug("start babel: %s %s", os.getcwd(), os.system('ls'))

    # prepare files for scwrl4
    babel_pdb_for_scwrl(wtpdb, proper_pdb=babeled_pdb, proper_fasta='wt.fasta')
    open(seq_s4_in, 'w').write(seq)

    # run scwrl4
    run_scwrl(pdb_s4_out, log_s4_out, babeled_pdb, seq_s4_in)

    # add deleted residues back into pdb
    post_process_scwrlpdb(resids, wtpdb, pdb_s4_out, new4q)

    # compute clash-score
    import clash as clash
    # clash_score, clash_resids = clash.clashscore_and_residues(resids, open(new4q).readlines())  # NOPEP8
    clash_score, clash_resids = clash.clashscore_and_residues(resids, new4q)  # NOPEP8

    if clash_score > 0:
        # dont work on residues that are marked as immutable:
        logger.debug('immutable: %s', immutable_resids)
        for immutable in immutable_resids:
            if immutable in clash_resids:
                clash_resids.remove(immutable)

        logger.info('Clash-Score was: %s, will now re-scrwl and allow to change residues %s', clash_score, sorted(clash_resids))

        def backup(fil):
            if os.path.basename(fil) != fil:
                raise ValueError('fil must be a filename, not a path', fil)
            shutil.move(fil, '.bkp'+fil)

        def restore(fil):
            if not os.path.isfile('.bkp'+fil):
                raise IOError('is not a file', '.bkp'+fil)
            shutil.move('.bkp'+fil, fil)

        def remove_backup(fil):
            if not os.path.isfile('.bkp'+fil):
                raise IOError('is not a file', '.bkp'+fil)
            os.remove('.bkp'+fil)

        # backing files up
        backup(log_s4_out)
        backup(pdb_s4_out)
        backup(seq_s4_in)
        backup(new4q)

        # list of residues that should be removed
        kill_residues = []

        # add clashing residues to sequence
        lseq = list(seq)
        for m in clash_resids:
            logger.debug("Residue clashes: %s %s %s", m, m-1, len(seq))
            if (m-1) >= len(lseq):
                # residue that comes after the protein sequence!
                # if solvent molecule, we should kill this residue in new4q.
                if m in immutable_resids:
                    logger.debug('keep residue number %s: is immutable', m)
                else:
                    kill_residues.append(m)
                continue
            lseq[m-1] = lseq[m-1].upper()
        newseq = ''.join(lseq)
        
        logger.debug('oldseq: %s', seq)
        logger.debug('newseq: %s', newseq)

        # overwrite sequence
        open(seq_s4_in, 'w').write(newseq)

        # re-run scrwl
        run_scwrl(pdb_s4_out, log_s4_out, babeled_pdb, seq_s4_in)
        post_process_scwrlpdb(clash_resids, wtpdb, pdb_s4_out, new4q)

        for resi in kill_residues:
            logger.debug('remove residue with number %s', resi)
            kill_residue(new4q, resi)

        # new_clash_score, new_clash_resids = clash.clashscore_and_residues(clash_resids, open(new4q).readlines())  # NOPEP8
        new_clash_score, new_clash_resids = clash.clashscore_and_residues(clash_resids, new4q)  # NOPEP8


        if new_clash_score < clash_score:
            logger.info('Clash-Score new: %s ==>  Keep.', new_clash_score)
            clash_resids = new_clash_resids
            clash_score = new_clash_score
            remove_backup(log_s4_out)
            remove_backup(pdb_s4_out)
            remove_backup(seq_s4_in)
            remove_backup(new4q)
        else:
            logger.info('Clash-Score new: %s ==>  Rollback!', new_clash_score)
            # restore backups
            restore(log_s4_out)
            restore(pdb_s4_out)
            restore(seq_s4_in)
            restore(new4q)
    else:
        logger.info('No clash detected')

    # use qprep5 to obtain smooth newq
    import qprep5
    newq, top = qprep5.create_top(qprep5inp, os.getcwd(), new4q)

    from fep import create_fep
    create_fep(wtpdb, wtfep, newq, outfep='mutant.fep')

    logger.debug('timing: %s s', round(time.time()-start, 2))


if __name__ == "__main__":
    import sys
    def usage():
        print(sys.argv[0], 'wt.pdb', 'wt.fep', 'wt.qpinp', 'sequence', ' [immutable resids]')
        print('or')
        print(sys.argv[0], 'wt.pdb', 'wt.fep', 'wt.qpinp', 'RESID:NEWAA', ' [immutable resids]')
        sys.exit(1)

    if len(sys.argv) < 5:
        usage()

    immutable_resids = []
    for resid in sys.argv[6:]:
        immutable_resids.append(int(resid))

    if ":" in sys.argv[4]:
        scwrl(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], immutable_resids)
    else:
        scwrl2(sys.argv[1], sys.argv[4], sys.argv[2], sys.argv[3], immutable_resids)


