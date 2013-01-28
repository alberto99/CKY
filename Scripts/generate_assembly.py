#! /usr/bin/env python
import assemble


st_1 = ['h1h2.rep.c1','h1h2.rep.c3','h1h2.rep.c4']
st_2 = ['h3.rep.c1','h3.rep.c2','h3.rep.c4']

i = 0
for s1 in st_1:
    for s2 in st_2:
        assemble.make_assembly(first_sse=s1,sec_sse=s2,loop_aa='DGLDEDG',output_name='template_%i.pdb' % i)
        i += 1

