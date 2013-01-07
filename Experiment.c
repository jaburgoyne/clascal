// CLASCAL (Experiment.c)
//
// Copyright (c) 2010 by John Ashley Burgoyne and the Royal Institute for the 
// Advancement of Learning (McGill University). All rights reserved.
//
// This source is adapted from Suzanne Winsberg's CLASCAL, version 7.01 (May
// 1993), written in FORTRAN 77. 
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice, 
//      this list of conditions, and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions, and the following disclaimer in the 
//      documentation and/or other materials provided with the distribution.
//
//   3. Neither the name of McGill University nor the names of its contributors
//      may be used to endorse or promote products derived from this software 
//      without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>
#include "inlines.h"
#include "StimulusSet.h"
#include "_StimulusSet.h"
#include "SubjectSet.h"
#include "_SubjectSet.h"
#include "Experiment.h"
#include "_Experiment.h"

struct _Experiment {
        char * restrict description;
        StimulusSet * restrict stimulusSet;
        SubjectSet * restrict subjectSet;
        double * restrict dissimilarities;
        size_t dataSize;
        double * squaredDistances;
};



/* N.B. Basic initialisation must occur before calling. */
static double * NewSquaredDistances(const Experiment * restrict self)
{
        if (!self->dissimilarities) return NULL;
        const size_t stimulusPairCount = StimulusPairCount(self->stimulusSet);
        const size_t subjectPairCount = SubjectPairCount(self->subjectSet);
        if (!subjectPairCount) return NULL;
        const SubjectPair * subjectPairs = SubjectPairs(self->subjectSet);
        double * squaredDistances;
        squaredDistances = SafeMalloc(subjectPairCount, sizeof(double));
        double * accumulator = SafeMalloc(stimulusPairCount, sizeof(double));
        for (size_t l = 0; l < subjectPairCount; l++) {
                cblas_dcopy((int)stimulusPairCount, 
                            (self->dissimilarities 
                             + stimulusPairCount * subjectPairs[l].i2), 
                            1, 
                            accumulator, 
                            1);
                cblas_daxpy((int)stimulusPairCount, 
                            -1.0, 
                            (self->dissimilarities 
                             + stimulusPairCount * subjectPairs[l].i1), 
                            1, 
                            accumulator, 
                            1);
                squaredDistances[l] = cblas_ddot((int)stimulusPairCount, 
                                                 accumulator, 
                                                 1, 
                                                 accumulator, 
                                                 1);
        }
        free(accumulator);
        return squaredDistances;
}

Experiment * NewExperiment(char * description,
                           StimulusSet * restrict stimulusSet,
                           SubjectSet * restrict subjectSet,
                           double * restrict dissimilarities)
{
        if (!stimulusSet || !subjectSet || !dissimilarities) return NULL;
        const size_t dataSize = SizeProduct(StimulusPairCount(stimulusSet),
                                            SubjectCount(subjectSet));
        Experiment * restrict self;
        if ((self = malloc(sizeof(Experiment)))) {
                self->description = description;
                self->stimulusSet = stimulusSet;
                self->subjectSet = subjectSet;
                self->dissimilarities = dissimilarities;
                self->dataSize = dataSize;
                self->squaredDistances = NewSquaredDistances(self);
        }
        return self;
}

void DeleteExperiment(Experiment * restrict self)
{
        if (self) {
                FreeAndClear(self->squaredDistances);
                FreeAndClear(self->dissimilarities);
                FreeAndClear(self->subjectSet);
                FreeAndClear(self->stimulusSet);
                FreeAndClear(self->description);
                free(self);
        }
}

const char * ExperimentDescription(const Experiment * self)
{
        return self ? self->description : NULL;
}

const StimulusSet * ExperimentStimulusSet(const Experiment * restrict self)
{
        return self ? self->stimulusSet : NULL;
}

const SubjectSet * ExperimentSubjectSet(const Experiment * restrict self)
{
        return self ? self->subjectSet : NULL;
}

const double * SubjectDissimilarities(const Experiment * restrict self)
{
        return self ? self->dissimilarities : NULL; 
}

size_t DataSize(const Experiment * self)
{
        return self ? self->dataSize : 0;
}

double * ExperimentSquaredDistances(const Experiment * self)
{
        return self ? self->squaredDistances : NULL;
}
