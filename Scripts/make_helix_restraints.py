#!/usr/bin/env python
# encoding: utf-8


h1 = [4, 18]

for start, end in [h1]:
    end = end - 4 + 1

    for i in range(start, end):
        print '%d \t O \t %d \t N \t 0.0 \t 3.5 \t 0.6 \t 0.0' % (i, i+4)

