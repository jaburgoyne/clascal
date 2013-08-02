// CLASCAL (ClassAssignment.c)
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
#include "Model.h"
#include "ClassAssignment.h"
#include "_ClassAssignment.h"

struct _ClassAssignment {
        const SubjectSet * restrict subjectSet;
        const Model * restrict model;
        double * restrict conditionalDistributions;
        size_t distributionsSize;
        double * restrict classSizes;
        double * restrict unconditionalDistribution;
};

static double * NewClassSizes(ClassAssignment * restrict self)
{
        const size_t classCount = ClassCount(self->model);
        if (!classCount) return NULL;
        const size_t subjectCount = SubjectCount(self->subjectSet);
        if (!subjectCount) return NULL;
        double * restrict classSizes = SafeCalloc(classCount, sizeof(double));
        for (size_t i = 0; i < subjectCount; i++)
                cblas_daxpy((int)classCount, 
                            1.0, 
                            self->conditionalDistributions + classCount * i, 
                            1, 
                            classSizes, 
                            1);
        return classSizes;
}

/* 
 * N.B.: This function may only be called after class sizes have been 
 * initialised. 
 */
static double * NewUnconditionalDist(const ClassAssignment * restrict self)
{
        if (!self->classSizes) return NULL;
        const size_t classCount = ClassCount(self->model);
        if (!classCount) return NULL;
        const size_t subjectCount = SubjectCount(self->subjectSet);
        if (!subjectCount) return NULL;
        const double scalingFactor = 1.0 / (double)subjectCount;
        double * restrict unconditionalDistribution;
        unconditionalDistribution = SafeMalloc(classCount, sizeof(double));
        cblas_dcopy((int)classCount, 
                    self->classSizes, 
                    1, 
                    unconditionalDistribution, 
                    1);
        cblas_dscal((int)classCount, 
                    scalingFactor, 
                    unconditionalDistribution, 
                    1);
        return unconditionalDistribution;
}

static double * 
NewNormalisedDistributions(const ClassAssignment * restrict self,
                           const double * restrict distributions)
{
        const size_t subjectCount = SubjectCount(self->subjectSet);
        const size_t classCount = ClassCount(self->model);
        double * restrict normalisedDistributions;
        normalisedDistributions = SafeMalloc(self->distributionsSize,
                                             sizeof(double));
        cblas_dcopy((int)self->distributionsSize, 
                    distributions, 
                    1, 
                    normalisedDistributions, 
                    1);
        for (size_t i = 0; i < subjectCount; i++) {
                // N.B.: All values must be non-negative.
                const double sum = cblas_dasum((int)classCount, 
                                               (normalisedDistributions
                                                + classCount * i), 
                                               1);
                cblas_dscal((int)classCount, 
                            1.0 / sum, 
                            normalisedDistributions + classCount * i, 
                            1);
        }
        return normalisedDistributions;
}

ClassAssignment * NewClassAssignment(const SubjectSet * restrict subjectSet,
                                     const Model * restrict model,
                                     const double * restrict distributions)
{
        if (!subjectSet || !model) return NULL;
        const size_t distributionsSize = SizeProduct(SubjectCount(subjectSet),
                                                     ClassCount(model));
        if (!IsConvertibleToInt(distributionsSize))
                ExitWithError("Too many classes to assign subjects");
        ClassAssignment * restrict self;
        if ((self = malloc(sizeof(ClassAssignment)))) {
                self->subjectSet = subjectSet;
                self->model = model;
                self->distributionsSize = distributionsSize;
                if (distributions) {
                        double * restrict normDists;
                        normDists = NewNormalisedDistributions(self,
                                                               distributions);
                        self->conditionalDistributions = normDists;
                } else {
                        double * restrict evenDists;
                        evenDists = SafeMalloc(distributionsSize, 
                                               sizeof(double));
                        for (size_t i = 0; i < distributionsSize; i++)
                                evenDists[i] = 1.0;
                        double * restrict normDists;
                        normDists = NewNormalisedDistributions(self, evenDists);
                        self->conditionalDistributions = normDists;
                        FreeAndClear(evenDists);
                }
                self->classSizes = NewClassSizes(self);
                self->unconditionalDistribution = NewUnconditionalDist(self);
        }
        return self;
}

void DeleteClassAssignment(ClassAssignment * restrict self)
{
        if (self) {
                FreeAndClear(self->unconditionalDistribution);
                FreeAndClear(self->classSizes);
                FreeAndClear(self->conditionalDistributions);
                FreeAndClear(self);
        }
}

const SubjectSet * 
ClassAssignmentSubjectSet(const ClassAssignment * restrict self)
{
        return self ? self->subjectSet : NULL;
}

const Model * ClassAssignmentModel(const ClassAssignment * restrict self)
{
        return self ? self->model : NULL;
}

const double * SubjectClassDistributions(const ClassAssignment * restrict self)
{
        return self ? self->conditionalDistributions : NULL;
}

const double * ClassSizes(const ClassAssignment * restrict self)
{
        return self ? self->classSizes : NULL;
}

const double * 
UnconditionalClassDistribution(const ClassAssignment * restrict self)
{
        return self ? self->unconditionalDistribution : NULL;
}

size_t DistributionsSize(const ClassAssignment * restrict self)
{
        return self ? self->distributionsSize : 0;
}
