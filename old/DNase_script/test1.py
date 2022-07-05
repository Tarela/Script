#!/usr/bin/env python
import numpy

a = range(400)
b = numpy.array(range(400))#numpy.array(a)
#c= sum(b)
#c=[]
for k in range(int(1e4)):
#    c.append(k)
   # for i in a:
   #     c=0
   #     c = b[i]

    [b[i] for i in a]


