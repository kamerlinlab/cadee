#!/usr/bin/env python

"""Shell Tool Launcher
Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""

import sys
import os


def main(args):
    def tool_usage(badtool=False):
        print()
        if badtool:
            print('Error')
            print('    Unknown Tool:', badtool)
        print('    Available Tools: lossy_repack(lr) | repair_simpack(rs)')
        sys.exit(1)

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
    else:
        tool_usage(badtool=subcmd)

    shell_command = shell_command + " " + " ".join(sys.argv)

    os.system(shell_command)
