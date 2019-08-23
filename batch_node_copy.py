from optparse import OptionParser
import subprocess as sub

from collections import defaultdict
import threading

import sys,os
from mpi4py import MPI
from sys import stderr

from  time import sleep
import difflib
import subprocess as sub

TAG_READY=1
TAG_DIE=2
TAG_DONE_JOB=3
TAG_DONE_COPY=4
TAG_RUN_JOB=5
TAG_ERROR=6


def mkdir(dir,file):
    ls_dir = os.listdir(dir)
    if(not(file_exists(ls_dir,file))):
        command = "mkdir %s/%s"%(dir,file)
        os.system(command)
    return "%s/%s"%(dir,file)


def unlink_all(to_unlink):
    for f in to_unlink:
        os.unlink(f)


def client_do_command(rank, proc_name, cmd, dry_run):
    if dry_run:
        print "CLIENT %s: rank %d, command:\n%s" % (proc_name, myrank, cmd)
        return 0
    else:
        p = sub.Popen(cmd,stdout=sub.PIPE,stderr=sub.PIPE,shell=True)
        output,errors = p.communicate()

        print "CLIENT %s: rank %d, o: %s e:%s\n%s"%(proc_name,myrank,output,errors,cmd)
        ret_code = p.returncode
        if ret_code != 0: print "".join(["********************\n" for i in xrange(1)])

        return ret_code


def do_rsync(myrank, proc_name, source, dest, pre_sync_commands, post_sync_commands, rsync_options, dry_run):
    if rsync_options is None:
        rsync_options = "-arzv"

    # Run the precopy comands then the rsync, then the post copy
    if pre_sync_commands != None:
        print "CLIENT %s: rank %d doing pre-sync jobs..."%(proc_name,myrank)
        client_do_command(myrank,proc_name,"%s"%(pre_sync_commands), dry_run)

    print "CLIENT %s: rank %d doing rsync..."%(proc_name,myrank)
    client_do_command(myrank, proc_name, "rsync %s %s %s" % (rsync_options, source, dest), dry_run)

    if post_sync_commands != None:
        print "CLIENT %s: rank %d doing post-sync jobs..."%(proc_name,myrank)
        client_do_command(myrank,proc_name,"%s"%(post_sync_commands), dry_run)


def client_loop(comm, myrank, proc_name, n_procs, source, dest, pre_sync_commands, post_sync_commands, rsync_options, dry_run):
    print "CLIENT %s: rank %d initialized..."%(proc_name,myrank)

    # Tell the boss I'm A-OK.
    comm.send({"rank":myrank,"proc_name":proc_name},dest=0, tag=TAG_READY)
    done_running = False

    while not done_running:
        status = MPI.Status()
        data = comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
        print "CLIENT %s: rank %d received a message - data: %s"%(proc_name,myrank,data)

        if status.tag==TAG_DIE:
            print "CLIENT %s: rank %d finishing..."%(proc_name,myrank)
            done_running=True
        elif status.tag==TAG_RUN_JOB:
            do_rsync(myrank, proc_name, source, dest, pre_sync_commands, post_sync_commands, rsync_options, dry_run)

            #tell the boss I'm done w/ my jobs
            comm.send(dest=0,tag=TAG_DONE_JOB)
        else:
            print "CLIENT %s: %d received unknown tag %d"%(proc_name,myrank,status.tag)

    sys.exit(0)


def host_loop(comm, myrank, proc_name, n_procs, max_nodes, source, dest, pre_sync_commands, post_sync_commands, rsync_options, dry_run):
    """
    Init
    """
    l_ranks_working = []
    l_ranks_done_working = []

    l_ranks_ready = []
    nodes_to_ranks=defaultdict(list)
    ranks_to_nodes={}

    print "initting nodes..."
    status=MPI.Status()
    for i in xrange(1,n_procs):
        data = comm.recv(source = i,tag=MPI.ANY_TAG,status=status)
        print "HOST %s: %d ready, %s, with data: %s"%(proc_name,i,status,data)
        rank=i
        node=data["proc_name"]
        l_ranks_ready.append(i)
        nodes_to_ranks[node].append(i)
        ranks_to_nodes[i] = node

    print "HOST: all ranks initialized"

    # Track unique nodes and make a copy for initial list of waiting nodes.
    unique_nodes = nodes_to_ranks.keys()
    waiting_nodes = unique_nodes[:]
    completed_nodes = []

    # Start jobs on max nodes or if there are fewer waiting nodes than the
    # maximum, start on all nodes.
    num_nodes_to_start = min(len(waiting_nodes), max_nodes)
    print "HOST: sending work to %i nodes" % num_nodes_to_start

    # Now for each unique node, tell one rank on that node to run the copy process.
    for i in xrange(num_nodes_to_start):
        node = waiting_nodes.pop()
        rank = nodes_to_ranks[node][0]
        comm.send("", tag=TAG_RUN_JOB, dest=rank)

    # Wait until all unique nodes are completed.
    while len(completed_nodes) != len(unique_nodes):
        status = MPI.Status()
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        if status.tag == TAG_DONE_JOB:
            print "HOST %s: rank %d received a message - JOB DONE from %d" % (proc_name, myrank, status.source)
            completed_nodes.append(ranks_to_nodes[status.source])
            print "HOST: %i nodes remain in queue" % len(waiting_nodes)

            # Send work to the next waiting node if any remains.
            if len(waiting_nodes) > 0:
                node = waiting_nodes.pop()
                rank = nodes_to_ranks[node][0]
                comm.send("", tag=TAG_RUN_JOB, dest=rank)

    print "HOST %s: rank %d ALL WORKER THREADS DONE"%(proc_name,myrank)
    print "HOST: cleaning up workers..."
    for i in xrange(1, n_procs):
        comm.send('',tag=TAG_DIE,dest=i)

    print "HOST: syncing data to host node"
    do_rsync(myrank, proc_name, source, dest, pre_sync_commands, post_sync_commands, rsync_options, dry_run)

    print "HOST: DONE"

    MPI.Finalize()
    sys.exit(0)


if __name__=="__main__":
    opts = OptionParser()
    opts.add_option('', '--source', dest='source', default=None)
    opts.add_option('', '--dest', dest='dest', default=None)
    opts.add_option('', '--pre_sync_commands', dest='pre_sync_commands', default=None)
    opts.add_option('', '--post_sync_commands', dest='post_sync_commands', default=None)
    opts.add_option('', '--max_nodes', dest='max_nodes', type="int", default=4,
                    help="limit the number of nodes to which files are simultaneously copied")
    opts.add_option('', '--rsync_options', dest='rsync_options', default=None,
                    help="override default rsync options")
    opts.add_option('', '--dry-run', action="store_true", dest='dry_run', default=False,
                    help="print commands to run on each node without actually running them")

    (o, args) = opts.parse_args()

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    n_procs = comm.Get_size()
    proc_name = MPI.Get_processor_name()

    if myrank==0:
        host_loop(comm, myrank, proc_name, n_procs, o.max_nodes, o.source, o.dest, o.pre_sync_commands, o.post_sync_commands, o.rsync_options, o.dry_run)
    else:
        client_loop(comm, myrank, proc_name, n_procs, o.source, o.dest, o.pre_sync_commands, o.post_sync_commands, o.rsync_options, o.dry_run)
