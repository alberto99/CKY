import os
import templates

def make_dir(name):
    os.mkdir(name)

def write_line(f,data):
    fo = open(f,'w')
    fo.write("%s\n" % data)
    fo.close()

def write_setup(f,seq,cycles):
    fo = open(f,'w')
    seq = str(seq)
    if seq[0] == seq[-1]:
        startup = templates.from_string
        extra_restraints = '\n'
    else:
        startup = templates.from_template
        extra_restraints = templates.extra_restraints
    txt = templates.setup_script.substitute(name='Seq%s' % seq,startup=startup,cycles=cycles,extra_restraints=extra_restraints)
    fo.write(txt)
    fo.close()

def make_master(name,seq):
    txt = ''


    seq = str(seq)
    if seq[0] == seq[-1]:
        if int(seq[0]) == 1:
            txt = '''#PBS -N CKY
#PBS -j oe
#PBS -A TG-MCB120052

#PBS -l walltime=24:00:00
#PBS -l nodes=10:ppn=3:gpus=3:shared

date
cd $PBS_O_WORKDIR

echo "nodefile="
cat $PBS_NODEFILE
echo "=end nodefile"

module load python
module load scipy
module load numpy
'''
    else:
        assemble = templates.make_assembly.format(name=name)
    txt += templates.master_script.format(name=name,assembly=assemble)
    return txt
