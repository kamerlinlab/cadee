#!/usr/bin/env python
from __future__ import print_function

"""Shell Tool Launcher
Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""

import sys
import os


def main(args, caller):
    def tool_usage(exitcode):
        print()
        print()
        print('Available Tools:')
        print()
        print('       repair_simpack (rs):')
        print('           Description: Unpack and Repack a simpack to fix errors.')
        print('           Caution:     WILL OVERWRITE ORIGINAL SIMPACK!')
        print('           Example:     {0} repair_simpack /full/path/to/simpacks/wt_0.tar'.format(caller))
        print()
        print('       lossy_repack (lr):')
        print('           Description: Utility to repack simpacks lossy, deleting dcd files.')
        print('                        Will apply to all simpacks in current working directory')
        print('           Example:     cd /full/path/to/simpacks; {0} lossy_repack'.format(caller))
        print('')
        print('')
        sys.exit(exitcode)

    if len(args) < 2:
        tool_usage()

    subcmd = args[1].lower()
    args.remove(subcmd)
    args.remove(args[0])

    shell_command = "/bin/bash "
    shell_command += os.path.dirname(os.path.abspath(__file__))

    if   subcmd == 'delete_tempfiles' or subcmd == 'dt':
        shell_command = os.path.join(shell_command, 'delete_tempfiles.sh')
    elif subcmd == 'lossy_repack'     or subcmd == 'lr':
        shell_command = os.path.join(shell_command, 'lossy_repack.sh')
    elif subcmd == 'repair_simpack'   or subcmd == 'rs':
        shell_command = os.path.join(shell_command, 'repair_simpack.sh')
    elif subcmd == '--help':
        tool_usage(0)
    else:
        print("Unknown command:", subcmd)
        tool_usage(1)

    shell_command = shell_command + " " + " ".join(sys.argv)

    os.system(shell_command)

if __name__ == "__main__":
    main(sys.argv, sys.argv[0])

