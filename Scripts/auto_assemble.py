#! /usr/bin/env python
import assemble
import shutil
import initial_setup
import sys
import os
import create_constraints
import glob

#Defaults
numb_cl = 3

def create_restraints():
    print os.getcwd()
    os.chdir('TEMPLATES')
    files = glob.glob('template_*.pdb')
    for f in files:
        name = '{}restraints'.format(f[0:-3])
        print name,f
        create_constraints.main(fi=f,fo=name)
    os.chdir('..')
    print os.getcwd()

def get_seq_ss(dirname):
    os.chdir('..')
    system = initial_setup.CreateSystem()
    system.make_ss_groups()
    os.chdir(dirname)
    return system

def generate_pairs(match,system):

    print system.seq
    print system.ss
    print system.sse

    os.chdir('TEMPLATES')
    #want to keep in order of importance: c0+c0 is more important than c0+c3. Have to do fair for all 
    #pairings that lead to this Fragment (123 --> 1+23 and 12+3 --> the c0+c0 has to be first for both
    k = 0
    for i in range(numb_cl):
        for j in range(numb_cl):
            for (f1,f2) in match:
                s1 = 'h%s.cluster%i' % (f1,i)
                s2 = 'h%s.cluster%i' % (f2,j)
                (tmp,end,ty) = system.sse[int(str(f1)[-1])-1]
                (start,tmp,ty) = system.sse[int(str(f2)[0])-1]
                #The numbering is amber like, starts at 1. For beggining of loop (end of first fragment)
                #we would do -1 and +1 to select the first "." instead of "H" or "E", so keep same number
                #For the end of loop, we should do -1 and -1; since the [ini:fin] is a [:) f(x), only do -1
                loop = str(system.seq)[end:start-1]
                print loop
                assemble.make_assembly(first_sse=s1,sec_sse=s2,loop_aa=loop,output_name='template_%i.pdb' % k)

                fo = ('template_%i.restraints' % k,'w')
                create_constraints.main(fi=s1,fo=fo)
                create_constraints.main(fi=s2,fo=fo)
                k += 1
    os.chdir('..')

def copy_files(dirname):
    name = dirname[9:]
    print name
    fragment = len(name)
    print fragment
    match = []
    for i in range(fragment-1):
        first = name[0:i+1]
        second = name[i+1:]
        for j in range(numb_cl):
            shutil.copy2('../Fragment_%s/h%s.cluster%i' % (first,first,j),'TEMPLATES/.')
            shutil.copy2('../Fragment_%s/h%s.cluster%i' % (second,second,j),'TEMPLATES/.')
        match.append( (first,second) )
    return (match)


def submit_queue(dirname):
    txt = '''#!/bin/bash
#$ -S /bin/bash
#$ -R yes
#$ -cwd
#$ -V
#$ -N Assemble
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
        match = copy_files(dirname)
        system = get_seq_ss(dirname)
        generate_pairs(match,system)
        #create_restraints()
    else:
        os.chdir(dirname)
        submit_queue(dirname)



if __name__ == '__main__':
    main()

