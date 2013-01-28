#! /usr/bin/env python
from zam import protein,geometry
import numpy
import sys

def make_assembly(first_sse='',sec_sse='',loop_aa='',output_name='template.pdb'):
    h1 = protein.Protein(first_sse)
    h2 = protein.Protein(sec_sse)
    h1.RemoveTermini()
    h2.RemoveTermini()
    ext = protein.Protein(Seq=loop_aa)
    
    #keep end position of the first sse and the length of loop to know 
    #range of loop to do MC on
    #want to take both ends, overlapping 1 Amino Acid between first chain and loop
    end = len(h1.Res)
    loop_l = len(ext.Res) - 1
    
    #Get the first residue in loop region
    ext_ini = ext.Res[0].StartAtom
    ext_end = ext.Res[0].StopAtom
    P2 = ext.Pos[ext_ini:ext_end,:]
    #Get the last residue in first helix
    h1_ini = h1.Res[-1].StartAtom
    h1_end = h1.Res[-1].StopAtom
    P1 = h1.Pos[h1_ini:h1_end,:]
    print len(P1)
    print len(P2)
    h1.WritePdb('h1.pdb')
    ext.WritePdb('ext.pdb')
    rotmat = geometry.RotMatRMSD(P1,P2,Center=False)
    
    ext.Pos = numpy.dot(ext.Pos,rotmat)
    
    #Now put the first residue of helix2 into the last loop residue
    h2_ini = h2.Res[0].StartAtom
    h2_end = h2.Res[0].StopAtom
    P1 = h2.Pos[h2_ini:h2_end,:]
    #loop
    ext_ini = ext.Res[-1].StartAtom
    ext_end = ext.Res[-1].StopAtom
    P2 = ext.Pos[ext_ini:ext_end,:]
    rotmat = geometry.RotMatRMSD(P2,P1,Center=False)
    
    h2.Pos = numpy.dot(h2.Pos,rotmat)
    
    ext.DelRes(0)
    #This fails with one residue in loop, so make a try/pass
    try:
        ext.DelRes(len(ext.Res)-1)
    except:
        #Have to take out the first residue on the next fragment if only one residue in loop
        #otherwise, repeat of 1 amino acid
        h2.DelRes(0)
    h1.Join(ext)
    h1.Join(h2)
    
    #Arrange backbone in loop
    #This fails with one residue in loop, so make a try/pass
    try:
        h1.Optimize(ResInd = range(end,end+loop_l) , N = 100, Ti = 0.001, Tf = 0.0001,FracPhiPsiMove = 1., FracChainMove = 0.5, BBScale = 0.1)
    except:
        pass
    #Arrange sidechains
    h1.Optimize(ResInd = None, N = 100, Ti = 0.001, Tf = 0.0001,FracPhiPsiMove = 0.01, FracChainMove = 0., BBScale = 0.1)
    h1.WritePdb(output_name)


def main():
    s1 = sys.argv[1]
    s2 = sys.argv[2]
    loop = sys.argv[3]
    out = sys.argv[4]
    make_assembly(first_sse=s1,sec_sse=s2,loop_aa=loop,output_name=out)
    

if __name__ == '__main__':
    main()
