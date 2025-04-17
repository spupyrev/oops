#!/usr/bin/env python3

import multiprocessing
import os, sys
import random
import time
import psutil
import warnings

from queue import Empty, Queue
from subprocess import PIPE, Popen
from threading import Thread
from time import sleep

warnings.filterwarnings("ignore", category=RuntimeWarning)

################################################################################
SLEEP = 5
CONSOLE_WIDTH = 140
NUM_THREADS = min(max(3 * multiprocessing.cpu_count() // 4, 1), 10)
PARALLEL_RUN = True

NUM_THREADS = 6
#PARALLEL_RUN = False

N=10

CMD="./oops -verbose=0 -i=/home/spupyrev/research/one_planar/data/cub18.g6 -Ccross2 -Cnic -part=%(thread_id)d/%(num_threads)d"
# CMD="./oops -verbose=0 -i=/home/spupyrev/research/one_planar/data/cub22-gir6.g6 -Ccross2 -part=%(thread_id)d/%(num_threads)d"
################################################################################


def enqueue_output(out, queue):
    for line in iter(out.readline, b""):
        queue.put(line)
    out.close()


def enqueue_outputs(proc, queue):
    # enqueue_output(proc.stdout, queue)
    enqueue_output(proc.stderr, queue)


def LOG(thread_id, thread_desc, msg):
    msg = msg.ljust(CONSOLE_WIDTH)
    sys.stdout.write("\x1b7\x1b[%dA%s%s\x1b8" % (NUM_THREADS - thread_id + 2, thread_desc, msg))
    sys.stdout.flush()


def start_thread(thread_id, thread_desc, cmd):
    log_file = "log_" + str(thread_id)
    with open(log_file, "w") as f:
        f.write("stderr for thread " + str(thread_id) + "\n")

    LOG(thread_id, thread_desc, "starting execution")
    sleep(random.random() * SLEEP)

    p = Popen(cmd, stdout=PIPE, stderr=PIPE, bufsize=1, shell=True)
    q = Queue()
    t = Thread(target=enqueue_outputs, args=(p, q))
    t.daemon = True  # thread dies with the program
    t.start()

    try:
        # read line without blocking
        while p.poll() is None:
            sleep(SLEEP)

            try:
                line = q.get_nowait()
            except Empty:
                line = ""
            else:  # got line
                with open(log_file, "a") as f:
                    f.write(line.decode())
                    f.write("\n")
                qs = 1
                last_line = line.strip()
                while not q.empty():
                    line = q.get().strip()
                    if len(line) > 0:
                        last_line = line
                    # logging
                    with open(log_file, "a") as f:
                        f.write(line.decode())
                        f.write("\n")
                    qs += 1

                last_line = last_line.decode().strip()
                LOG(thread_id, thread_desc, last_line)                  

        if p.returncode == 0:
            with open(log_file) as f:
                for line in f:
                    pass
                last_line = line
            LOG(thread_id, thread_desc, "\033[92mok\033[0m {}".format(last_line.strip('\n')))
            return 0
        else:
            if p.returncode == -2:
                # ctrl+C
                LOG(thread_id, thread_desc, '\033[91m%s\033[0m' % ('failed with exit code ' + str(p.returncode)))
            else:
                # LOG(thread_id, thread_desc, '\033[91m%s\033[0m' % ('failed with exit code ' + str(p.returncode)))
                pass
    except Exception as e:
        LOG(thread_id, thread_desc, "\033[91m%s\033[0m" % ("process failed with error message '" + str(e) + "'"))
    return 1


def main():
    # remove previous logs
    for thread_id in range(0, 100):
        log_file = "log_" + str(thread_id)
        if os.path.exists(log_file):
            os.remove(log_file)

    # print system stats
    print("system info:")
    print("  cpu count: {}".format(psutil.cpu_count()))
    print("  cpu free : {:2.0f} %".format(100 - 100 * psutil.getloadavg()[1] / os.cpu_count()))
    print("  RAM avail: {:2.0f} GB".format(psutil.virtual_memory().available / (1024 ** 3)))
    print("  RAM total: {:2.0f} GB".format(psutil.virtual_memory().total / (1024 ** 3)))

    print("using %d (%s) threads for command '%s'" % (NUM_THREADS, "parallel" if PARALLEL_RUN else "sequential", CMD))
    print("=" * CONSOLE_WIDTH)
    for thread_id in range(0, NUM_THREADS):
        print("initializing thread {}...".format(thread_id))
    print("=" * CONSOLE_WIDTH)

    random.seed(123)
    seeds = []
    for thread_id in range(0, NUM_THREADS):
        seeds.append(random.randint(0, 1000000000))

    start_time = time.time()
    threads = []
    for thread_id in range(0, NUM_THREADS):
        seed = seeds[thread_id]

        # prepare the description
        thread_desc = "Thread {:2d}".format(thread_id + 1)
        if "(thread_id)" in CMD:
            if NUM_THREADS > 10:
                thread_desc += " (thread_id={:2d}; num_threads={})".format(thread_id, NUM_THREADS)
            else:
                thread_desc += " (thread_id={}; num_threads={})".format(thread_id, NUM_THREADS)
        if "(seed)" in CMD:
            thread_desc += " (seed={})".format(seed)
        thread_desc = thread_desc.ljust(30)
        thread_desc = "\033[90m" + thread_desc + ": \033[0m"

        # prepare the command
        cmd = (CMD) % {
            "seed": seed,
            "thread_id": thread_id,
            "num_threads": NUM_THREADS,
            "n": N,
        }

        # start the thread
        t = Thread(
            target=start_thread,
            args=(thread_id + 1, thread_desc, cmd),
        )
        t.daemon = True  # thread dies with the program
        t.start()
        threads.append(t)
        if not PARALLEL_RUN:
            t.join()

    for t in threads:
        t.join()

    end_time = time.time()
    print("computation completed in {:.0f} sec".format(end_time - start_time))


################################################################################
################################ MAIN ##########################################
################################################################################
if __name__ == "__main__":
    main()
