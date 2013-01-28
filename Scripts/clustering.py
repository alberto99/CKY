#! /usr/bin/env python
#Given a directory should do:
#1. gathering of the trajectory
#2. get topology
#3. create ptraj script and cluster
#4. Process ptraj output and slect structures 1,2,3 for next stage, give name according to fragment and order.

import shutil
import sys
import os

#Defaults
keep_clusters = 3

def process_clusters(dirname,keep_clusters=keep_clusters):
    name = dirname[9:]
    fi = open('h%s.txt' % name)
    ok = 0
    cluster = []
    for l in fi:
        if 'Clustering: divide' in l:
            ok = 1
            continue
        if not ok:
            continue
        #Cluster    0: has   249 points, occurence 0.144
        if 'occurence' in l:
            c = l.split()
            cluster.append((c[1][0:-1],c[3]))
        else:
            break
    #sort clusters by population
    best = sorted(cluster, key=lambda cl: cl[1],reverse=True)
    for i in range(keep_clusters):
        (c,p) = best[i]
        shutil.copy2('h%s.rep.c%s' % (name,c),'h%s.cluster%i' % (name,i))

def create_ptraj_script(dirname):
    name = dirname[9:]
    txt = '''ptraj topol.top<<EOF >& err_ptraj
    trajin traj.crd 200 100000 1
    rmsd first mass  @CA
    cluster out h%s representative pdb \ 
       average pdb means clusters 5 rms @CA
EOF
    ''' % name
    fo = open('ptraj.sh','w')
    fo.write(txt)
    fo.close()
    os.system('sh ptraj.sh')
    os.unlink('ptraj.sh')
    os.unlink('err_ptraj')

def pdb_to_traj():
    txt = '''set mol [mol new traj.pdb type pdb waitfor all]
    set sel [atomselect top "all"]
    animate write crd traj.crd sel $sel
    quit'''
    fo = open('vmd.tcl','w')
    fo.write(txt)
    fo.close()
    os.system('vmd  -dispdev text -eofexit < vmd.tcl')
    os.unlink('vmd.tcl')

def submit_queue(dirname):
    txt = '''#!/bin/bash
#$ -S /bin/bash
#$ -R yes
#$ -cwd
#$ -V
#$ -N Clust
#$ -j y
#$ -q cpu.q

module add epd

export PATH=/vault/tools/vmd/bin:$PATH

export LD_LIBRARY_PATH=/vault/tools/vmd/lib/vmd:$LD_LIBRARY_PATH

#disable using GPU
export VMDCUDANODISPLAYGPUS=1
export VMDNOCUDA=1 


%s QUEUE %s
''' % (sys.argv[0],sys.argv[1])
    fo = open('clust.sh','w')
    fo.write(txt)
    fo.close()
    os.system('qsub clust.sh')


def main():
    dirname = sys.argv[1]
    if dirname == 'QUEUE':
        dirname = sys.argv[2]
        os.system('pwd')
        #os.system('~/src/springs/bin/gather.py 0 0')
        os.system('~/src/springs/bin/topol.py')
        pdb_to_traj()
        create_ptraj_script(dirname)
        process_clusters(dirname)
    else:
        os.chdir(dirname)
        submit_queue(dirname)



if __name__ == '__main__':
    main()

