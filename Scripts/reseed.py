#! /usr/bin/env python
#Given a directory should do:
#1. gathering of the trajectory
#2. get topology
#3. create ptraj script and cluster
#4. Process ptraj output and slect structures 1,2,3 for next stage, give name according to fragment and order.

import shutil
import sys
import os
import numpy
import glob
import create_constraints
import templates
from zam import protein

#Defaults
keep_clusters = 3

def create_setup_file():
    fo = open('setup.py','w')
    txt = templates.reseed_template
    fo.write(txt)
    fo.close()

def create_restraints():
    os.chdir('TEMPLATES')
    files = glob.glob('template_*.pdb')
    for f in files:
        name = '{}restraints'.format(f[0:-3])
        create_constraints.generate_contacts(protein.Protein(f),fo=open(name,'w'),offset=0,wide=4,fce=0.1)
    os.chdir('..')

def create_backup(dirname):
    os.mkdir(dirname)
    files = ['traj.pdb.gz','traj.pdb.gz.backup','results.hd5','results.hd5.backup','restart.dat','restart.dat.backup','run_remd.log','TEMPLATES']
    for f in files:
        shutil.move(f,dirname)
    shutil.copy('setup.py',dirname)
    os.mkdir('TEMPLATES')

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

def create_pdb():
    txt = '''ptraj topol.top<<EOF >& err_ptraj
    trajin traj.crd 
    trajout TEMPLATES/template_ pdb
EOF
    '''
    fo = open('ptraj.sh','w')
    fo.write(txt)
    fo.close()
    os.system('sh ptraj.sh')
    os.unlink('ptraj.sh')
    os.unlink('err_ptraj')
    os.chdir('TEMPLATES')
    for i,f in enumerate(glob.glob('template_*')):
        shutil.move(f,'template_%s.pdb' % i)
    os.chdir('..')


def create_ptraj_script(dirname,trajin='trajin traj.crd 200 100000 1\n'):
    name = dirname[9:]
    txt = '''ptraj topol.top<<EOF >& err_ptraj
    {trajin}
    rmsd first mass  @CA
    cluster out h{name} representative pdb \ 
       average pdb means clusters 5 rms @CA
EOF
    '''.format(trajin=trajin,name=name)
    fo = open('ptraj.sh','w')
    fo.write(txt)
    fo.close()
    os.system('sh ptraj.sh')
    os.unlink('ptraj.sh')
    os.unlink('err_ptraj')

def pdb_to_traj(number):
    txt = '''set mol [mol new traj.pdb type pdb waitfor all]
    set sel [atomselect top "all"]
    set num_steps [molinfo top get numframes]
    set last [expr $num_steps-%s]
    animate delete beg 0 end $last top
    animate write crd traj.crd sel $sel
    quit''' % (number)
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

def submit_keeneland(frame,backdir):
    txt='''#PBS -N clust
#PBS -j oe
#PBS -A TG-MCB120052

#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1

date
cd $PBS_O_WORKDIR

echo "nodefile="
cat $PBS_NODEFILE
echo "=end nodefile"

module load python
module load scipy
module load numpy

export PYTHONPATH=~/src/springs:$PYTHONPATH
export PYTHONPATH=~/src/Python:$PYTHONPATH

module add epd
#disable using GPU
export VMDCUDANODISPLAYGPUS=1
export VMDNOCUDA=1 


%s QUEUE %s %s
''' % (sys.argv[0],frame,backdir)
    fo = open('clust.sh','w')
    fo.write(txt)
    fo.close()
    os.system('qsub clust.sh')

def analyze_goap(directory,top=1000,criteria='dfire'):
    path = os.getcwd()
    os.chdir(directory)
    files = glob.glob('pdb.*')
    files = '\n'.join(files)
    txt = "/nics/d/home/alberto3/src/goap-alone\n"
    txt += files
    fo = open('goap.in','w')
    fo.write(txt)
    fo.close()
    os.system('/nics/d/home/alberto3/src/goap-alone/goap < goap.in > results.txt')
    return select_top(top=top,criteria=criteria,path=path)

def select_top(top=1000,path='..',criteria='dfire'):
    data = numpy.loadtxt('results.txt', dtype=({'names':['i','name','both','dfire','goap'],'formats':[numpy.int, 'S15', numpy.float,numpy.float,numpy.float]}))

    a = numpy.argsort(data['dfire'])[0:1000]
    b = numpy.argsort(data['goap'])[0:1000]
    c = numpy.argsort(data['both'])[0:1000]
    alla =numpy.concatenate((a,b,c))
    a = numpy.unique(alla)
    c = ['trajin PDB/%s' % b for b in data['name'][a]]
    os.chdir(path)
    txt = '\n'.join(c)
    return txt


def main():
    #Num_Frames is the number of structures we are going to keep for reseeding
    num_frames = sys.argv[1]
    #Back_dir is the name of the directory were we are going to keep the old data
    back_dir = sys.argv[2]
    if num_frames == 'QUEUE':
        num_frames = sys.argv[2]
        back_dir = sys.argv[3]
        os.system('pwd')
        os.system('zcat traj.pdb.gz > traj.pdb')
        #os.system('~/src/springs/bin/gather.py 0 0')
        os.system('~/src/springs/bin/topol.py')
        pdb_to_traj(num_frames)
        #Create Backup folder for done trajectory
        create_backup(back_dir)
        #create pdbs and restraints
        create_pdb()
        create_restraints()
        #seed = int(back_dir.split("_")[-1])+1
        #create_setup_file(seed=seed)
        create_setup_file()
    else:
        #submit_queue(dirname)
        submit_keeneland(num_frames,back_dir)



if __name__ == '__main__':
    main()

