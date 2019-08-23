#!/usr/bin/env python
"""
Run rsync as an MPI job.
"""
import logging
from mpi4py import MPI
import optparse
import os
import subprocess
import sys
import time

# Setup logging.
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Log to stdout.
sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
sh.setFormatter(formatter)
logger.addHandler(sh)

# Define MPI constants.
MASTER=0
TAG_WORK=11
TAG_DONE=12
TAG_DIE=13


def get_next_work_item(command, options, input_file, parameters):
    if options.command_list:
        work = input_file
    else:
        parameters["dummy_in"] = input_file
        work = command % parameters

    return work


def start_mpi(command, options=None):
    """
    Run a given command through MPI. Input and output parameters are optionally
    defined in ``options``, a optparse options instance.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    node_name = MPI.Get_processor_name()
    number_of_nodes = comm.Get_size()
    nodes = xrange(1, number_of_nodes)
    status = MPI.Status()

    if rank == MASTER:
        logger.info("Master node is %s", node_name)

        if options.dry_run:
            logger.info("Dry run mode: no commands will be executed on nodes")

        logger.info("Start MPI for command: %s", command)
        logger.info("Got %s nodes", number_of_nodes)

        # Master nodes send out commands.
        if options.input_file:
            input_files = [line.strip()
                           for line in open(options.input_file, "r")]
            logger.debug("Found input files: %s", input_files)
        elif options.input_dir:
            input_files = [os.path.join(options.input_dir, filename)
                           for filename in os.listdir(options.input_dir)]
        else:
            input_files = []

        # Prepare parameters.
        parameters = {}
        if hasattr(options, "output_dir"):
            parameters["dummy_out"] = options.output_dir

        # Send work to nodes.
        working_nodes = []
        for node in nodes:
            if len(input_files) > 0:
                # Send work to each node as long as there is any left.
                work = get_next_work_item(command, options, input_files.pop(0), parameters)
                comm.send(work, dest=node, tag=TAG_WORK)
                working_nodes.append(node)
            else:
                # Kill any remaining nodes when there isn't any work left.
                comm.send("Die", dest=node, tag=TAG_DIE)

        # Send remaining work to nodes as they finish.
        while input_files:
            # Wait for any node to respond.
            result = comm.recv(
                source=MPI.ANY_SOURCE,
                tag=MPI.ANY_TAG,
                status=status
            )

            # Send next input to the node that finished.
            work = get_next_work_item(command, options, input_files.pop(0), parameters)

            comm.send(work, dest=status.source, tag=TAG_WORK)

        # Wait for working nodes to finish now that all work has been sent.
        for node in working_nodes:
            result = comm.recv(
                source=MPI.ANY_SOURCE,
                tag=MPI.ANY_TAG,
                status=status
            )
            logger.info("Node %s is done", status.source)

        # Shutdown all working nodes.
        for node in working_nodes:
            comm.send(None, dest=node, tag=TAG_DIE)

        logger.info("Master node is done.")
    else:
        # Non-master nodes process commands.
        while True:
            command = comm.recv(
                source=MASTER,
                tag=MPI.ANY_TAG,
                status=status
            )

            if status.tag == TAG_DIE:
                logger.info("Rank %s (node %s) complete", rank, node_name)
                break
            else:
                logger.info("Rank %s (node %s) running command: %s", rank, node_name, command)
                # Execute command on the filesystem if this isn't a dry run.
                if options.dry_run:
                    return_value = 0
                    run_time = 0
                else:
                    start_time = time.time()
                    return_value = subprocess.call(command, shell=True)
                    run_time = time.time() - start_time

                if return_value != 0:
                    logger.error(
                        "Command on rank %s (node %s) returned non-zero exit code: %s",
                        rank,
                        node_name,
                        return_value
                    )
                else:
                    logger.info("Rank %s (node %s) finished command (%s) in %.2f seconds", rank, node_name, command, run_time)

                comm.send("Done", dest=MASTER, tag=TAG_DONE)


if __name__ == "__main__":
    usage = """%prog [options] command"""
    parser = optparse.OptionParser()
    parser.add_option("--input_file", dest="input_file",
                      help="path to a text file with a list of commands or file paths to process on each line")
    parser.add_option("--input_dir", dest="input_dir",
                      help="path to a directory for which each file will be processed as input")
    parser.add_option("--output_dir", dest="output_dir",
                      help="path to a directory where output should be placed")
    parser.add_option("--suffix", dest="suffix",
                      help="suffix to append to each output file")
    parser.add_option("--commands",
                      action="store_true", dest="command_list", default=False,
                      help="input file contains one complete command per line")
    parser.add_option("--dry-run",
                      action="store_true", dest="dry_run", default=False,
                      help="Process all MPI functions without running any commands on the nodes.")
    (options, args) = parser.parse_args()

    if len(args) < 1 and not options.command_list:
        parser.error("You must specify a command.")
    else:
        if options.command_list:
            command = None
        else:
            command = args[0]

        start_mpi(command, options)
