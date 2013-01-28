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
    else:
        startup = templates.from_template
    txt = templates.setup_script.substitute(name='Seq%s' % seq,startup=startup,cycles=cycles)
    fo.write(txt)
    fo.close()

def make_master(name,seq):
    txt = ''
    seq = str(seq)
    if seq[0] == seq[-1]:
        pass
    else:
        txt += templates.make_assembly.format(name=name)
    txt += templates.master_script.format(name=name)
    return txt
