#! /usr/bin/env python
import numpy
from matplotlib import pyplot

aa = numpy.loadtxt('replica_trace.dat')


col = ['k-','g-','y-','b-','r-','m-','c-','k-','g-','y-','b-','r-','m-','c-','k-','g-','y-','b-','r-','m-','c-','k-','g-','y-','b-','r-','m-','c-','k-','g-','y-','b-','r-','m-','c-']
pyplot.figure()
ax = pyplot.subplot(111)
x = range(0,len(aa[0,:]))

for i in range(0,len(aa[:,0])):
    pyplot.plot(x,aa[i,:],col[i])

pyplot.savefig('remd.pdf',format='pdf')
