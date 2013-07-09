import string

from_template = 'ladder.init_from_templates( glob.glob(\'TEMPLATES/template_*.pdb\') )'
from_string = '''sequence = rest_parse.get_fasta_sequence( open('sequence.dat').read() )
        ladder.init_from_sequence(sequence)'''

extra_restraints = '''
        contact_files = glob.glob('TEMPLATES/template*.restraints')
        rests = []
        for c in contact_files:
            try:
                rests.append(restraints.RestraintCombiner(rest_parse.get_contact_restraints( open(c).read() ) ))
            except:
                pass
        if len(rests) > 0:
            ladder.add_restraints(rests,
            restraints.BinaryLowestECollection, accuracy=1./len(rests),
            force_scaler=contact_scaler)
        '''

setup_script = '''#!/usr/bin/env python

import numpy
import glob
import rest_parse
import system
import restraints


def add_restraint(start_index, delta_i, delta_j, add_helix):
    extended_bounds = numpy.load('bounds_E_5.npy')
    helical_bounds = numpy.load('bounds_H_5.npy')

    i = start_index + delta_i
    j = start_index + delta_j

    if add_helix:
        min_dist = helical_bounds[delta_i, delta_j, 0]
        max_dist = helical_bounds[delta_i, delta_j, 1]
    else:
        min_dist = extended_bounds[delta_i, delta_j, 0]
        max_dist = extended_bounds[delta_i, delta_j, 1]

    rest = restraints.DistanceRestraint(
            i, 'CA',
            j, 'CA',
            min_dist,
            max_dist,
            k = 1.5,
            quadratic_range = 2.0 )
    print '    Added restraint {}, {}, [{}, {}]'.format(i, j, min_dist, max_dist)
    return rest


def get_secondary_restraints():
    ss_string = open('ss.dat').read().strip()
    n_res = len(ss_string)

    for i in range(n_res):
        assert ss_string[i] in ['H', 'E', '.']

    rests = []

    for start in range(0, n_res - 4):
        end = start + 5

        add_helix = False
        add_extended = False

        if ss_string[start:end].count('E') >= 4:
            print '{} adding extended from {} to {}.'.format(
                    ss_string[start:end], start, end-1 )
            add_extended = True

        elif ss_string[start:end].count('H') >= 4:
            print '{} adding helix from {} to {}.'.format(
                    ss_string[start:end], start, end-1 )
            add_helix = True

        if add_extended or add_helix:
            #rest1 = add_restraint(start+1, 0, 1, add_helix)
            #rest2 = add_restraint(start+1, 1, 2, add_helix)
            #rest3 = add_restraint(start+1, 2, 3, add_helix)
            #rest4 = add_restraint(start+1, 3, 4, add_helix)

            #rest5 = add_restraint(start+1, 0, 2, add_helix)
            #rest6 = add_restraint(start+1, 1, 3, add_helix)
            #rest7 = add_restraint(start+1, 2, 4, add_helix)

            rest8 = add_restraint(start+1, 0, 3, add_helix)
            rest9 = add_restraint(start+1, 1, 4, add_helix)

            rest10 = add_restraint(start+1, 0, 4, add_helix)

            rests.append( restraints.RestraintCombiner(
                #[rest1, rest2, rest3, rest4, rest5, rest6, rest7, rest8, rest9, rest10] ) )
                [rest8, rest9, rest10] ) )

    return rests

def sort_by_co(dist_rests):
    def get_co(dist):
        return abs(dist.i - dist.j)
    return sorted(dist_rests, key=get_co, reverse=True)

def sort_and_setup(dist_rests):
    sorted_dists = sort_by_co(dist_rests)
    n = len(sorted_dists)

    w = 0.25

    for (i, dist) in enumerate(sorted_dists):
        # first we map i into [0, 1.0 - w)
        alpha_min = float(i) / float(n - 1) * (1.0 - w)
        alpha_max = alpha_min + w
        print dist.i, dist.j, abs(dist.i - dist.j), alpha_min, alpha_max
        dist._force_scaler = restraints.NonLinearScaling(alpha_min, alpha_max, 3.5)

    return sorted_dists


ladder = system.Ladder()

with ladder.setup():
    with ladder.setup_parameters():
        ${startup}
        ladder.system_name = '${name}'
        ladder.cycles = ${cycles}
        ladder.ramp_steps = 20
        ladder.force_field = 'leaprc.ff12SB'
        ladder.t_min = 270
        ladder.t_max = 400
        ladder.backup_frequency = 100
        ladder.ps_per_stage = 10
        ladder.igb_model = 8
        ladder.use_amap = True
        ladder.n_replicas = 30

        ladder.alpha_t_min = 0.0
        ladder.alpha_t_max = 1.0

        ladder.adaptive_remd_active = True
        ladder.adaptive_remd_burn_in = 50
        ladder.adaptive_remd_length = 50
        ladder.adaptive_remd_decay = 2.0

    with ladder.setup_restraints():
        ss2 = rest_parse.get_secondary_restraints( open('ss.dat').read(), force_const=1.0 )
        ladder.add_restraints(ss2, restraints.ConstantCollection)

        contact_scaler = restraints.NonLinearScaling(0.0, 1.0, 8.0)

        #Fragments of 5 and lower for sec structure:
        ss_restraints2 = rest_parse.get_secondary_restraints( open('ss.dat').read() )
        ladder.add_restraints(ss_restraints2, restraints.BinaryLowestECollection, accuracy=0.7)

        ${extra_restraints}

        confinement_rest = restraints.get_confinement_restraints( len(ladder.sequence), 30.0, 1.0 )
        ladder.add_restraints(confinement_rest, restraints.ConstantCollection)
'''

setup_script = string.Template(setup_script)


make_assembly = '''~/src/CKY/Scripts/auto_assemble.py QUEUE {name}'''
master_script = '''
cd {name}
{assembly}
chmod +x setup.py
export PYTHONPATH=~/src/springs:$PYTHONPATH
./setup.py
~/src/springs/bin/run_remd.py
zcat traj.pdb.gz > traj.pdb
module add epd
export PATH=/vault/tools/vmd/bin:$PATH
export LD_LIBRARY_PATH=/vault/tools/vmd/lib/vmd:$LD_LIBRARY_PATH
#disable using GPU
export VMDCUDANODISPLAYGPUS=1
export VMDNOCUDA=1 
~/src/CKY/Scripts/clustering.py QUEUE {name}
cd ..
'''

reseed_template = '''#!/usr/bin/env python

import numpy
import glob
import rest_parse
import system
import restraints
import shutil
import random


def make_hydroph_groups():
        ss =  open('ss.dat').read()
        extended = 0
        helical = 0
        sse = []
        for i,l in enumerate(ss):
            if l not in "HE.":
                continue
            if l in 'E':
                if extended:
                    continue
                else:
                    start = i
                    extended = 1
            if l in 'H':
                if helical:
                    continue
                else:
                    start = i
                    helical = 1
            if l not in 'E' and extended:
                end = i
                sse.append((start+1,end+1))
                extended = 0
            if l not in 'H' and helical:
                end = i
                sse.append((start+1,end+1))
                helical = 0
        print sse
        return sse

def add_restraint(start_index, delta_i, delta_j, add_helix):
    extended_bounds = numpy.load('bounds_E_5.npy')
    helical_bounds = numpy.load('bounds_H_5.npy')

    i = start_index + delta_i
    j = start_index + delta_j

    if add_helix:
        min_dist = helical_bounds[delta_i, delta_j, 0]
        max_dist = helical_bounds[delta_i, delta_j, 1]
    else:
        min_dist = extended_bounds[delta_i, delta_j, 0]
        max_dist = extended_bounds[delta_i, delta_j, 1]

    rest = restraints.DistanceRestraint(
            i, 'CA',
            j, 'CA',
            min_dist,
            max_dist,
            k = 1.5,
            quadratic_range = 2.0 )
    print '    Added restraint {}, {}, [{}, {}]'.format(i, j, min_dist, max_dist)
    return rest


def get_secondary_restraints():
    ss_string = open('ss.dat').read().strip()
    n_res = len(ss_string)

    for i in range(n_res):
        assert ss_string[i] in ['H', 'E', '.']

    rests = []

    for start in range(0, n_res - 4):
        end = start + 5

        add_helix = False
        add_extended = False

        if ss_string[start:end].count('E') >= 4:
            print '{} adding extended from {} to {}.'.format(
                    ss_string[start:end], start, end-1 )
            add_extended = True

        elif ss_string[start:end].count('H') >= 4:
            print '{} adding helix from {} to {}.'.format(
                    ss_string[start:end], start, end-1 )
            add_helix = True

        if add_extended or add_helix:
            #rest1 = add_restraint(start+1, 0, 1, add_helix)
            #rest2 = add_restraint(start+1, 1, 2, add_helix)
            #rest3 = add_restraint(start+1, 2, 3, add_helix)
            #rest4 = add_restraint(start+1, 3, 4, add_helix)

            #rest5 = add_restraint(start+1, 0, 2, add_helix)
            #rest6 = add_restraint(start+1, 1, 3, add_helix)
            #rest7 = add_restraint(start+1, 2, 4, add_helix)

            rest8 = add_restraint(start+1, 0, 3, add_helix)
            rest9 = add_restraint(start+1, 1, 4, add_helix)

            rest10 = add_restraint(start+1, 0, 4, add_helix)

            rests.append( restraints.RestraintCombiner(
                #[rest1, rest2, rest3, rest4, rest5, rest6, rest7, rest8, rest9, rest10] ) )
                [rest8, rest9, rest10] ) )

    return rests

def sort_by_co(dist_rests):
    def get_co(dist):
        return abs(dist.i - dist.j)
    return sorted(dist_rests, key=get_co, reverse=True)

def sort_and_setup(dist_rests):
    sorted_dists = sort_by_co(dist_rests)
    n = len(sorted_dists)

    w = 0.25

    for (i, dist) in enumerate(sorted_dists):
        # first we map i into [0, 1.0 - w)
        alpha_min = float(i) / float(n - 1) * (1.0 - w)
        alpha_max = alpha_min + w
        print dist.i, dist.j, abs(dist.i - dist.j), alpha_min, alpha_max
        dist._force_scaler = restraints.NonLinearScaling(alpha_min, alpha_max, 3.5)

    return sorted_dists


ladder = system.Ladder()

with ladder.setup(cluster='keeneland'):
    with ladder.setup_parameters():
        ladder.init_from_templates( glob.glob('TEMPLATES/template_*.pdb') )
        ladder.system_name = 'Seq123'
        ladder.cycles = 5000
        ladder.ramp_steps = 20
        ladder.force_field = 'leaprc.ff12SB'
        ladder.t_min = 270
        ladder.t_max = 400
        ladder.backup_frequency = 100
        ladder.ps_per_stage = 10
        ladder.igb_model = 8
        ladder.use_amap = True
        ladder.n_replicas = 30

        ladder.alpha_t_min = 0.0
        ladder.alpha_t_max = 1.0

        ladder.adaptive_remd_active = True
        ladder.adaptive_remd_burn_in = 50
        ladder.adaptive_remd_length = 50
        ladder.adaptive_remd_decay = 2.0

    with ladder.setup_restraints():
        ss2 = rest_parse.get_secondary_restraints( open('ss.dat').read(), force_const=1.0 )
        ladder.add_restraints(ss2, restraints.ConstantCollection)

        contact_scaler = restraints.NonLinearScaling(0.0, 1.0, 8.0)

        #Fragments of 5 and lower for sec structure:
        ss_restraints2 = rest_parse.get_secondary_restraints( open('ss.dat').read() )
        ladder.add_restraints(ss_restraints2, restraints.BinaryLowestECollection, accuracy=0.5)

        destination = open('all_restraints.dat', 'wb')
        for filename in glob.glob('TEMPLATES/template*.restraints'):
            shutil.copyfileobj(open(filename, 'rb'), destination)
        destination.close()

        long_restraints = rest_parse.get_distance_bound_restraints( open('all_restraints.dat').read() )
        random.shuffle(long_restraints)
        n_long = float(len(long_restraints))
        if n_long > 500:
            long_restraints = long_restraints[:500]
            n_long = float(len(long_restraints))
        acc = 10.0/n_long
        ladder.add_restraints(
                long_restraints,
                restraints.BinaryLowestECollection, accuracy=acc,
                force_scaler=contact_scaler)

        sse = make_hydroph_groups()

        n_res = len(ladder.sequence)
        for i in range(len(sse)):
            for j in range(i+1,len(sse)):
                g1s,g1e = sse[i]
                g2s,g2e = sse[j]
                print g1s,g1e,g2s,g2e
                hydrophobic,l1,l2 = restraints.get_hydrophobic_contact_restraints(ladder.sequence, group_1=range(g1s,g1e), group_2=range(g2s,g2e), force_constant=0.1)
                ladder.add_restraints(hydrophobic, restraints.BinaryLowestECollection, accuracy=0.1,force_scaler=contact_scaler)


        confinement_rest = restraints.get_confinement_restraints( len(ladder.sequence), 30.0, 1.0 )
        ladder.add_restraints(confinement_rest, restraints.ConstantCollection)
'''

