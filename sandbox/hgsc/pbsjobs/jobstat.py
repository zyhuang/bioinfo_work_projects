import sys
from subprocess import PIPE, Popen
import time
import datetime


def jobstat(sleep_time=10):

    cmd = 'qstat -tu zhuoyih | grep zhuoyih'
    proc = Popen(cmd, stdout=PIPE, universal_newlines=True, shell=True)

    job_stat = {}
    for line in proc.stdout:

        # 1779539[3].sug-moab     zhuoyih     scavenge calc_pop_distanc  18108     1      1    5gb       --  R  00:00:00
        jobid, uname, qname, jname, foo, nnode, ncpu, mem, foo, status, rtime = line.rstrip().split()

        jobid = jobid.split('[')[0].split('.')[0]
        if jobid not in job_stat:
            job_stat[jobid] = {'T':0, 'C':0, 'R':0, 'Q':0, 'E':0, 'name': jname}
        job_stat[jobid][status] += 1
        job_stat[jobid]['T'] += 1

    print(datetime.datetime.today())
    for job,data in sorted(job_stat.items()):
        print('job: {}\tname: {}\ttotal: {}\tqueue: {}\trunning: {}\tcomplete: {}'
              .format(job, data['name'], data['T'], data['Q'], data['R'], data['C']))
    time.sleep(sleep_time)
    print()


if __name__ == '__main__':

    if len(sys.argv) == 1:
        sleep_time = 10
    else:
        sleep_time = int(sys.argv[1])

    while True:
        jobstat(sleep_time)
