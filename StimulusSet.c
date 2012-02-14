// CLASCAL (StimulusSet.c)
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
#include "inlines.h"
#include "StimulusSet.h"
#include "_StimulusSet.h"

struct _StimulusSet {
        size_t stimulusCount;
        char * * restrict stimulusNames;
        size_t pairCount;
        StimulusPair * restrict stimulusPairs;
};

/* N.B.: The stimulus pair must be otherwise initialised before calling. */
static StimulusPair * NewStimulusPairs(const StimulusSet * restrict self)
{
        if (!self) return NULL;
        StimulusPair * restrict stimulusPairs;
        stimulusPairs = SafeMalloc(self->pairCount, sizeof(StimulusPair));
        size_t m = 0;
        for (size_t k = 1; k < self->stimulusCount; k++) {
                for (size_t j = 0; j < k; j++) {
                        stimulusPairs[m].j = j;
                        stimulusPairs[m].k = k;
                        m++;
                }
        }
        return stimulusPairs;
}        

StimulusSet * NewStimulusSet(size_t stimulusCount,
                             char * * restrict stimulusNames)
{
        if (stimulusCount < 2 || !stimulusNames) return NULL;
        StimulusSet * restrict self;
        if ((self = malloc(sizeof(StimulusSet)))) {
                self->stimulusCount = stimulusCount;
                self->stimulusNames = stimulusNames;
                self->pairCount = (SizeProduct(stimulusCount, stimulusCount - 1) 
                                   / 2);
                self->stimulusPairs = NewStimulusPairs(self);
        }
        return self;
}

void DeleteStimulusSet(StimulusSet * restrict self)
{
        if (self) {
                for (size_t j = 0; j < self->stimulusCount; j++) {
                        FreeAndClear(self->stimulusNames[j]);
                }
                FreeAndClear(self->stimulusNames);
                FreeAndClear(self->stimulusPairs);
                free(self);
        }
}

size_t StimulusCount(const StimulusSet * restrict self)
{
        return self ? self->stimulusCount : 0;
}

const char * StimulusNameForIndex(const StimulusSet * restrict self, 
                                  size_t index) 
{
        return (self && index < self->stimulusCount
                ? self->stimulusNames[index]
                : NULL);
}

size_t StimulusPairCount(const StimulusSet * restrict self)
{
        return self ? self->pairCount : 0;
}

const StimulusPair * StimulusPairs(const StimulusSet * restrict self)
{
        return self ? self->stimulusPairs : NULL;
}

size_t PairNumberForStimulusPair(const StimulusSet * restrict self,
                                 const StimulusPair * restrict pair)
{
        if (!self || !pair) return 0;
        const size_t j = pair->j < pair->k ? pair->j : pair->k;
        const size_t k = pair->j < pair->k ? pair->k : pair->j;
        if (j >= self->stimulusCount || k >= self->stimulusCount || j == k)
                return self->pairCount;
        const size_t offset = SizeProduct(k, (k - 1)) / 2;
        return SizeSum(j, offset);
}
