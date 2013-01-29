#!/usr/bin/env python
# encoding: utf-8


import numpy
from zam import protein


excluded_residues = []



def generate_distances(p,fo='distances.dat'):
    rest = open(fo, 'w')

    n_res = len(p)

    for i in range(n_res):
        if (i+1) in excluded_residues:
            continue

        ind_i = i + 1
        print i, 'of', n_res

        if p.Seq[i] == 'GLY':
            name_i = 'CA'
        else:
            name_i = 'CB'

        for j in range(i+8, n_res):
            if (j+1) in excluded_residues:
                continue

            if p.Seq[j] == 'GLY':
                name_j = 'CA'
            else:
                name_j = 'CB'

            ind_j = j + 1

            pos_1 = p.Pos[p.AtomInd(AtomName=name_i, ResNum=i), :]
            pos_2 = p.Pos[p.AtomInd(AtomName=name_j, ResNum=j), :]

            d = numpy.linalg.norm(pos_1 - pos_2)

            if d < 10.0:

                print >>rest, '%d \t %s \t %d \t %s \t %f \t %f \t 0.0' % (
                        ind_i, name_i, ind_j, name_j, 10.0, 0.10 )


def generate_torsions(p):
    rest_10 = open('torsions10.dat', 'w')
    rest_30 = open('torsions30.dat', 'w')

    for i in range( len(p) ):
        if (i+1) in excluded_residues:
            continue

        ind_i = i + 1
        print i, 'of', len(p)

        phi, psi = p.PhiPsi(i)
        print >>rest_10, '%d \t %f \t 10.0 \t %f \t 10.0 \t 0.3 0.0' % (
                ind_i, phi, psi)
        print >>rest_30, '%d \t %f \t 30.0 \t %f \t 30.0 \t 0.3 0.0' % (
                ind_i, phi, psi)

def generate_ss(p):
    rest = open('secondary_box.dat', 'w')

    n_res = len(p)

    for i in range(1, n_res-1):
        if p.Seq[i] == 'GLY' or p.Seq[i] == 'PRO':
            continue
        res_i = i + 1
        print >>rest, '%d \t -105 45 140 30 0.3 0.0' % res_i
        print >>rest, '%d \t -80 30 -30 30 0.3 1.5' % res_i


def main(fi='template.pdb',fo='distances.dat'):
    p = protein.Protein(fi)

    generate_distances(p,fo)
    #generate_torsions(p)
    #generate_ss(p)

if __name__ == '__main__':
    main()
