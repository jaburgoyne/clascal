// CLASCAL (Model.c)
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
#include "Model.h"
#include "_Model.h"

struct _Model {
        size_t classCount;
        size_t dimensionCount;
        SpecificityType specificityType;
        size_t weightsSize;
        size_t extendedWeightsSize;
};

Model * NewModel(size_t dimensionCount,
                 size_t classCount,
                 SpecificityType specificityType)
{
        if (!classCount) return NULL;
        // TODO: This check can be removed when star trees are supported.
        if (!dimensionCount && specificityType) 
                ExitWithError("Star-tree models are not yet supported.");
        const size_t weightsSize = SizeProduct(dimensionCount, classCount);
        const size_t extendedWeightsSize = ((specificityType 
                                             == GlobalSpecificities)
                                            ? SizeSum(weightsSize, classCount)
                                            : weightsSize);        
        if (!IsConvertibleToInt(extendedWeightsSize))
                ExitWithError("Model is too large to process");
        Model * restrict self;
        if ((self = malloc(sizeof(Model)))) {
                self->dimensionCount = dimensionCount;
                self->classCount = classCount;
                self->specificityType = specificityType;
                self->weightsSize = weightsSize;
                self->extendedWeightsSize = extendedWeightsSize;
        }
        return self;
}

void DeleteModel(Model * restrict self)
{
        free(self);
}

size_t ClassCount(const Model * restrict self) 
{
        return self ? self->classCount : 0;
}

size_t DimensionCount(const Model * restrict self)
{
        return self ? self->dimensionCount : 0;
}

SpecificityType ModelSpecificityType(const Model * restrict self)
{
        return self ? self->specificityType : NoSpecificities;
}

size_t WeightsSize(const Model * restrict self)
{
        return self ? self->weightsSize : 0;
}

size_t ExtendedWeightsSize(const Model * restrict self)
{
        return self ? self->extendedWeightsSize : 0;
}

bool IsSpatial(const Model * restrict self)
{
        return self ? self->dimensionCount || self->specificityType : false;
}
