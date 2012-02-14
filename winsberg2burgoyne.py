#!/usr/bin/env python
#
#  winsberg2burgoyne.py
#  CLASCAL
#
# Copyright (c) 2010 by John Ashley Burgoyne and the Royal Institute for the 
# Advancement of Learning (McGill University). All rights reserved.

import sys

winsberg = open(sys.argv[1], 'r')
title = winsberg.readline().strip()
runParameters = winsberg.readline().split()
logicalConstants = winsberg.readline().split()
specificityConstants = winsberg.readline().split()
iterationConstants = winsberg.readline().split()
dataFormat = winsberg.readline().split()
data = winsberg.read().replace('D', 'E').split()
winsberg.close()

burgoyne = open(sys.argv[2], 'w')
burgoyne.write('{\n')
burgoyne.write('    "Description" : "' + title + '",\n\n')
subjectCount = int(runParameters[0])
stimulusCount = int(runParameters[1])
burgoyne.write('    "StimulusNames" : [\n')
for j in range(stimulusCount):
    burgoyne.write('        "Stimulus' + ("%02u" % j) + '"')
    if j < stimulusCount - 1: burgoyne.write(',')
    burgoyne.write('\n')
burgoyne.write('    ],\n\n')
burgoyne.write('    "SubjectData" : {\n\n')
for i in range(subjectCount):
    burgoyne.write('        "Subject' + ("%03u" % i) + '" : [\n')
    for k in range(1, stimulusCount):
        burgoyne.write('            [')
        for j in range(0, k):
            burgoyne.write("%.3f" % float(data.pop(0)))
            if j < k - 1 : burgoyne.write(', ')
        burgoyne.write(']')
        if k < stimulusCount - 1: burgoyne.write(',')
        burgoyne.write('\n')
    burgoyne.write('        ]')
    if i < subjectCount - 1: burgoyne.write(',')
    burgoyne.write('\n\n')
burgoyne.write('    }\n')
burgoyne.write('}\n')
burgoyne.close()
