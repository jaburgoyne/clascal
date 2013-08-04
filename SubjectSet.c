// CLASCAL (SubjectSet.c)
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
#include "SubjectSet.h"
#include "_SubjectSet.h"

struct _SubjectSet {
        size_t subjectCount;
        char * * restrict subjectNames;
        size_t pairCount;
        SubjectPair * restrict subjectPairs;
};

/* N.B. The subject pair must be otherwise initialised before calling. */
static SubjectPair * NewSubjectPairs(const SubjectSet * restrict self)
{
        if (!self) return NULL;
        SubjectPair * restrict subjectPairs;
        subjectPairs = SafeMalloc(self->pairCount, sizeof(SubjectPair));
        size_t l = 0;
        for (size_t i2 = 1; i2 < self->subjectCount; i2++) {
                for (size_t i1 = 0; i1 < i2; i1++) {
                        subjectPairs[l].i1 = i1;
                        subjectPairs[l].i2 = i2;
                        l++;
                }
        }
        return subjectPairs;
}

SubjectSet * NewSubjectSet(size_t subjectCount,
                           char * * restrict subjectNames)
{
        if (!subjectCount || !subjectNames) return NULL;
        if (!IsConvertibleToInt(subjectCount))
                ExitWithError("Too many subjects to process");
        SubjectSet * restrict self;
        if ((self = malloc(sizeof(SubjectSet)))) {
                self->subjectCount = subjectCount;
                self->subjectNames = subjectNames;
                self->pairCount = (SizeProduct(subjectCount, subjectCount - 1)
                                   / 2);
                self->subjectPairs = NewSubjectPairs(self);
        }
        return self;
}

void DeleteSubjectSet(SubjectSet * restrict self) {
        if (self) {
                FreeAndClear(self->subjectPairs);
                for (size_t i = 0; i < self->subjectCount; i++) {
                        FreeAndClear(self->subjectNames[i]);
                        
                }
                FreeAndClear(self->subjectNames);
                free(self);                
        }
}

size_t SubjectCount(const SubjectSet * restrict self)
{
        return self ? self->subjectCount : 0;
}

const char * SubjectNameForIndex(const SubjectSet * restrict self, 
                                 size_t index)
{
        return ((self && index < self->subjectCount)
                ? self->subjectNames[index]
                : NULL);
}

size_t SubjectPairCount(const SubjectSet * restrict self)
{
        return self ? self->pairCount : 0;
}

const SubjectPair * SubjectPairs(const SubjectSet * restrict self)
{
        return self ? self->subjectPairs : NULL;
}

size_t PairNumberForSubjectPair(const SubjectSet * restrict self,
                                const SubjectPair * restrict pair)
{
        if (!self) return 0;
        if (!pair) return self->pairCount;
        const size_t i1 = pair->i1 < pair->i2 ? pair->i1 : pair->i2;
        const size_t i2 = pair->i1 < pair->i2 ? pair->i2 : pair->i1;
        if (i1 >= self->subjectCount || i2 >= self->subjectCount || i1 == i2)
                return self->pairCount;
        const size_t offset = SizeProduct(i2, (i2 - 1)) / 2;
        return SizeSum(i1, offset);
}

