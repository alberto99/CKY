import string

from_template = 'ladder.init_from_templates( glob.glob(\'TEMPLATES/template_*.pdb\') )'
from_string = '''sequence = rest_parse.get_fasta_sequence( open('sequence.dat').read() )
        ladder.init_from_sequence(sequence)'''

extra_restraints = '''
        contact_files = glob.glob('TEMPLATES/template*.restraint')
        rests = []
        for c in contact_files:
            rests.append(restraints.RestrainCombiner(rest_parse.get_contact_restraints( open(c).read() ) ))
        ladder.add_restraints(rests,
        restraints.BinaryMonteCarloCollection, accuracy=1./len(rests),
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
        ladder.add_restraints(ss_restraints2, restraints.BinaryMonteCarloCollection, accuracy=0.7)

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

