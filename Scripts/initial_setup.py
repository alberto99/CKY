#! /usr/bin/env python
# encoding: utf-8
import copy
import os
import cky
import generate_files as gf

"""
look for an ss.dat and sequence.dat
given those two files:
    determine the number of secondary structure elements
    crease parse_tree directories
    create a setup.py, ss.dat, sequence.dat, and helix_dist.dat
"""

class CreateSystem(object):
    def __init__(self,ss='ss.dat',seq='sequence.dat'):
        self.seq = open(seq).read().strip()
        self.ss = open(ss).read().strip()
        self.sse = None

    def complete_pipeline(self):
        self.make_ss_groups()
        self.make_cky_tree()
        self.make_data()

    def make_cky_tree(self):
        self.levels = len(self.sse)
        self.cky = cky.SetupCKY(self.levels)

    def make_data(self):
        #write sequence and sse.dat files for every directory
        self.master =''
        for t in self.cky.states:
            self.get_sequence(t)
            self.get_cycles(t)
            dirname = 'Fragment_%s' % t
            gf.make_dir(dirname)
            gf.make_dir('%s/TEMPLATES' % dirname)
            gf.write_line('%s/ss.dat' % dirname,self.tmp_ss)
            gf.write_line('%s/sequence.dat' % dirname,self.tmp_seq)
            gf.write_setup('%s/setup.py' % dirname,t,self.cycles)
            self.create_helix_dist(dirname)
            self.master += gf.make_master(dirname,t)
        print self.master

    def get_cycles(self,name):
        poss_cycles = [500,1500,5000,5000,5000,5000]
        length = len(str(name))
        if length > 6:
            self.cycles=5000
        else:
            self.cycles=poss_cycles[length-1]

    def get_sequence(self,name):
        #A name is something like 123, makes reference t which sse are being used
        #self.sse has the numbering like amber: starts 1 to whichever.
        #However, when geting seq and ss we want to get the numbering -1. 
        #the range [ini:end] returns [ini:end); what we want is to extend in one position
        #up and down (get overlapping loop) and make the ss = "." for this element. So what
        #we want is [ini-1:end+1]
        name = str(name)
        self.tmp_sse = []
        for char in name:
            self.tmp_sse.append(self.sse[int(char)-1])
        if name[0] == name[-1]:
            #have only one sse.
            ini,end,typ = self.sse[int(name)-1]
        else:
            ini,temp,typ = self.sse[int(name[0])-1]
            temp,end,typ = self.sse[int(name[-1])-1]
        #Get the end elements in a 0 based array instead of 1 based
        ini = ini - 1
        end = end - 1
        #Correct for ends
        if  str(1) in str(name):
            ini = 0
        else:
            ini = ini-1
        if str(self.levels) in str(name): 
            end = len(self.ss)
        else:
            end = end+1

        self.tmp_ss = self.ss[ini:end+1]
        if ini != 0:
            self.tmp_ss = ".%s" % self.tmp_ss[1:]
        if end != len(self.ss):
            self.tmp_ss = "%s." % self.tmp_ss[:-1]
        self.tmp_seq = self.seq[ini:end+1]
        self.tmp_offset = ini-1



    def create_helix_dist(self,dirname):
        #Create 1-4 distance restraints for helical sse
        fe = open('%s/extended.dat' % dirname,'w')
        fh = open('%s/helix_dist.dat' % dirname,'w')
        for start, end,ss_type in self.tmp_sse:
            start = start - self.tmp_offset
            end = end - self.tmp_offset
            if ss_type == 'E':
                #for i in range(start,end-2+1)
                    #at 4.5 sd of 1.5, range 4-6 no energy; exclude too much of beta
                    #print '%d \t O \t %d \t N \t 4.0 \t 6.0 \t 0.6 \t 0.0' % (i, i+2)
                for i in range(start,end-3+1):
                    #at 8.0 sd of 2.0 under, 1 over
                    fe.write('%d \t O \t %d \t N \t 6.2 \t 9.0 \t 0.6 \t 0.0\n' % (i, i+3))
                for i in range(start,end-4+1):
                    #at 11.0 sd of 1.5 under, 1 over
                    fe.write('%d \t O \t %d \t N \t 8.2 \t 13.0 \t 0.6 \t 0.0\n' % (i, i+4))
            if ss_type == 'H':
                for i in range(start, end-4+1):
                    fh.write('%d \t O \t %d \t N \t 0.0 \t 4.0 \t 0.6 \t 0.0\n' % (i, i+4))
        fe.close()
        fh.close()

    def make_ss_groups(self):
            extended = 0
            helical = 0
            sse = []
            for i,l in enumerate(self.ss):
                if l not in "HE.":
                    continue
                if l not in 'E' and extended:
                    end = i
                    sse.append((start+1,end,'E'))
                    extended = 0
                if l not in 'H' and helical:
                    end = i
                    sse.append((start+1,end,'H'))
                    helical = 0
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
            print sse
            self.sse = sse


def main():
    system = CreateSystem()
    system.complete_pipeline()



if __name__ == '__main__':
    main()
