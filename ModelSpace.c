// CLASCAL (ModelSpace.c)
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
#include <string.h>
#include <Accelerate/Accelerate.h>
#include "inlines.h"
#include "StimulusSet.h"
#include "_StimulusSet.h"
#include "Model.h"
#include "_Model.h"
#include "ModelSpace.h"
#include "_ModelSpace.h"

// The following macros are only to make computing degrees of freedom
// somewhat more readable.

#define SP SizeProduct
#define SS SizeSum
#define SD(x, y) _SizeDifference(x, y, __FILE__, __LINE__, __func__)
static inline size_t _SizeDifference(size_t x, 
                                     size_t y, 
                                     const char * fileName,
                                     int lineNumber,
                                     const char * functionName)
{
        if (x < y)
                _ExitWithError("Difference of sizes would be negative",
                               fileName,
                               lineNumber,
                               functionName);
        return x - y;
}

struct _ModelSpace {
        const StimulusSet * restrict stimulusSet;
        const Model * restrict model;
        size_t coordinatesSize;
        size_t specificitiesSize;
        size_t differencesSize;
        size_t distancesSize;
        size_t parameterCount;
        double * restrict weights;
        double * restrict coordinates;
        double * restrict specificities;
        double * restrict pairSpecificities; // rows classes & cols stim pairs
        double * restrict differences; // rows stimulus pairs & cols dimensions
        double * restrict squaredDifferences;
        double * restrict distances; // rows classes & cols stimulus pairs
};

/* N.B.: Basic initialisation must be complete before calling. */
static size_t ParameterCountValue(const ModelSpace * restrict self)
{
        const size_t J = StimulusCount(self->stimulusSet);
        const size_t T = ClassCount(self->model);
        const size_t R = DimensionCount(self->model);
        const SpecificityType specificityType = ModelSpecificityType(self
                                                                     ->model);
        if (!IsSpatial(self->model)) return SS(SP(T,  SP(J, SD(J, 1)) / 2), T);
        // No break statements are necessary due to return statements.
        switch(specificityType) {
                case NoSpecificities:
                        return ((T == 1)
                                ? SD(SS(T, SP(J, R)), SP(R, SS(R, 1)) / 2)
                                : SD(SS(SS(T, SP(J, R)), SP(R, T)), SP(R, 2)));
                case GlobalSpecificities:
                        return ((T == 1)
                                ? SD(SS(SS(T, SP(J, R)), J),
                                     SP(R, SS(R, 1)) / 2) 
                                : SD(SS(SS(SS(T, SP(J, R)), SP(SS(R, 1), T)), 
                                        J),
                                     SS(SP(R, 2), 1))); 
                case ClassSpecificities:
                        return ((T == 1)
                                ? SD(SS(SS(T, SP(J, R)), J),
                                     SP(R, SS(R, 1)) / 2)
                                : SD(SS(SS(SS(T, SP(J, R)), SP(R, T)), 
                                        SP(J, T)),
                                     SP(R, 2)));
                default:
                        return 0;
        }
}

/* N.B.: Basic initialisation must be complete before calling. */
static double * NewUnfoldedSpecificities(const ModelSpace * restrict self)
{
        if (!self->specificities) return NULL;
        const size_t stimulusCount = StimulusCount(self->stimulusSet);
        if (!stimulusCount) return NULL;
        const size_t pairCount = StimulusPairCount(self->stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(self->stimulusSet);
        if (!stimulusPairs) return NULL;                
        const size_t classCount = ClassCount(self->model);
        if (!classCount) return NULL;
        const SpecificityType specType = ModelSpecificityType(self->model);
        const size_t unfoldedRows = (specType == GlobalSpecificities
                                     ? 1
                                     : classCount);
        const size_t unfoldedSize = SizeProduct(unfoldedRows, pairCount);
        double * restrict pairSpecificities  = SafeMalloc(unfoldedSize, 
                                                          sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                cblas_dcopy((int)unfoldedRows, 
                            self->specificities + stimulusPairs[m].j, 
                            (int)stimulusCount, 
                            pairSpecificities + m, 
                            (int)pairCount);
                cblas_daxpy((int)unfoldedRows, 
                            1.0, 
                            self->specificities + stimulusPairs[m].k, 
                            (int)stimulusCount, 
                            pairSpecificities + m, 
                            (int)pairCount);
        }
        return pairSpecificities;
}

/* N.B.: Basic initialisation must be complete before calling. */
static double * NewDifferences(const ModelSpace * restrict self)
{
        if (!self->coordinates) return NULL;
        const size_t stimulusCount = StimulusCount(self->stimulusSet);
        if (!stimulusCount) return NULL;
        const size_t pairCount = StimulusPairCount(self->stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(self->stimulusSet);
        if (!stimulusPairs) return NULL;                
        const size_t dimensionCount = DimensionCount(self->model);
        if (!dimensionCount) return NULL;
        double * restrict differences = SafeMalloc(self->differencesSize, 
                                                   sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                cblas_dcopy((int)dimensionCount, 
                            (self->coordinates 
                             + dimensionCount * stimulusPairs[m].k), 
                            1, 
                            differences + dimensionCount * m,
                            1);
                cblas_daxpy((int)dimensionCount,
                            -1.0, 
                            (self->coordinates 
                             + dimensionCount * stimulusPairs[m].j),
                            1,
                            differences + dimensionCount * m, 
                            1);
        }
        return differences;
}

/*
 * N.B.: Basic initialisation and differences must be complete before calling.
 */
static double * NewSquaredDifferences(const ModelSpace * restrict self)
{
        if (!self->differences) return NULL;
        double * restrict squaredDifferences = SafeMalloc(self->differencesSize, 
                                                          sizeof(double));
        cblas_dcopy((int)self->differencesSize,
                    self->differences,
                    1,
                    squaredDifferences,
                    1);
        cblas_dtbmv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)self->differencesSize, 
                    0,
                    self->differences, 
                    1, 
                    squaredDifferences,
                    1); // trick to effect elementwise multiplication     
        return squaredDifferences;
}

/* 
 * N.B.: Basic initialisation, differences, and unfolded specificities must be 
 * complete before calling. 
 */
static double * NewDistances(const ModelSpace * restrict self)
{
        if (!self->squaredDifferences) return NULL;
        const size_t stimulusCount = StimulusCount(self->stimulusSet);
        if (!stimulusCount) return NULL;
        const size_t pairCount = StimulusPairCount(self->stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(self->stimulusSet);
        if (!stimulusPairs) return NULL;                
        const size_t dimCount = DimensionCount(self->model);
        if (!dimCount) return NULL;
        const size_t classCount = ClassCount(self->model);
        if (!classCount) return NULL;
        if (!self->weights && classCount > 1) return NULL;
        const SpecificityType specificityType = ModelSpecificityType(self->
                                                                     model);
        if (specificityType != NoSpecificities && !self->pairSpecificities)
                return NULL;
        double * restrict distances = SafeCalloc(self->distancesSize, 
                                                 sizeof(double));
        // Compute the weighted sums of the squared differences.
#ifdef WINSBERG_LEGACY
        // Winsberg has a check for weights equal to zero only in the case
        // where there are no specificities. 
        if (!specificityType && classCount > 1) {
                for (size_t t = 0; t < classCount; t++) {
                        for (size_t m = 0; m < pairCount; m++) {
                                for (size_t r = 0; r < dimCount; r++) {
                                        const double w 
                                        = self->weights[classCount * r + t];
                                        const double d2
                                        = self->squaredDifferences[dimCount * m
                                                                   + r];
                                        distances[pairCount * t +m] 
                                        += (w == 0.0 ? 0.00001 * d2 : w * d2);
                                }
                        }
                }
                for (size_t i = 0; i < self->distancesSize; i++)
                        distances[i] = sqrt(distances[i]);
                return distances;                
        }
#endif
        if (classCount > 1) {
                cblas_dgemm(CblasRowMajor, 
                            CblasTrans, 
                            CblasTrans, 
                            (int)classCount, 
                            (int)pairCount, 
                            (int)dimCount, 
                            1.0, 
                            self->weights, 
                            (int)classCount, 
                            self->squaredDifferences, 
                            (int)dimCount, 
                            0.0, 
                            distances, 
                            (int)pairCount);
        } else {
                // The for loop seems better than allocating a vector of ones 
                // and calling a matrix multiply. Because the values have been 
                // squared, the absolute value is not a problem.
                for (size_t m = 0; m < pairCount; m++)
                        distances[m] = cblas_dasum((int)dimCount, 
                                                   (self->squaredDifferences
                                                    + dimCount * m), 
                                                   1);
        }
        // Add specificities (if necessary).
        switch (specificityType) {
                case ClassSpecificities:
                        cblas_daxpy((int)self->distancesSize, 
                                    1.0, 
                                    self->pairSpecificities, 
                                    1, 
                                    distances, 
                                    1);
                        break;
                case GlobalSpecificities:
                        if (classCount > 1)
                                cblas_dger(CblasRowMajor, 
                                           (int)classCount, 
                                           (int)pairCount, 
                                           1.0, 
                                           (self->weights 
                                            + WeightsSize(self->model)), 
                                           1, 
                                           self->pairSpecificities, 
                                           1, 
                                           distances, 
                                           (int)pairCount);
                        else 
                                cblas_daxpy((int)pairCount, 
                                            1.0, 
                                            self->pairSpecificities, 
                                            1, 
                                            distances, 
                                            1);                
                        break;
                default:
                        break;
        }
        // Take the square root of all elements.
        for (size_t i = 0; i < self->distancesSize; i++)
                distances[i] = sqrt(distances[i]);
        return distances;
}

/* 
 * Returns the average value over classes for each weight. 
 */
static double * NewWeightMeans(const ModelSpace * restrict self,
                               const double * weights)
{
        if (!weights) return NULL;
        const size_t classCount = ClassCount(self->model);
        const size_t weightCount = ((ModelSpecificityType(self->model) 
                                     == GlobalSpecificities)
                                    ? SizeSum(DimensionCount(self->model), 1)
                                    : DimensionCount(self->model));        
        double * weightSums = SafeCalloc(weightCount, sizeof(double));
        for (size_t t = 0; t < classCount; t++)
                cblas_daxpy((int)weightCount, 
                            1.0 / (double)classCount, 
                            weights + t, 
                            (int)classCount, 
                            weightSums, 
                            1);
        return weightSums;
}

/* N.B.: Basic initialisation must occur before calling. */
static double * NewNormalisedWeights(const ModelSpace * restrict self,
                                     const double * restrict weightMeans, 
                                     const double * restrict weights)
{
        if (!weights) return NULL;
        const size_t extendedWeightsSize = ExtendedWeightsSize(self->model);
        double * restrict normalisedWeights = SafeMalloc(extendedWeightsSize, 
                                                         sizeof(double));
        cblas_dcopy((int)extendedWeightsSize, 
                    weights, 
                    1, 
                    normalisedWeights, 
                    1);
        if (weightMeans) {
                const size_t classCount = ClassCount(self->model);
                const size_t weightCount = ((ModelSpecificityType(self->model)
                                             == GlobalSpecificities)
                                            ? SizeSum(DimensionCount(self
                                                                     ->model),
                                                      1)
                                            : DimensionCount(self->model));
                for (size_t r = 0; r < weightCount; r++) 
                        cblas_dscal((int)classCount, 
                                    1.0 / weightMeans[r], 
                                    normalisedWeights + classCount * r,
                                    1); 
        }
        return normalisedWeights;
}

/* N.B.: Basic initialisation must occur before calling. */
static double * NewNormalisedCoordinates(const ModelSpace * restrict self,
                                         const double * restrict weightMeans,
                                         const double * restrict coordinates)
{
        if (!coordinates) return NULL;
        const size_t stimulusCount = StimulusCount(self->stimulusSet);
        const size_t dimensionCount = DimensionCount(self->model);
        double * restrict normalisedCoordinates;
        normalisedCoordinates = SafeMalloc(self->coordinatesSize,
                                           sizeof(double));
        cblas_dcopy((int)self->coordinatesSize, 
                    coordinates, 
                    1, 
                    normalisedCoordinates, 
                    1);
        double * accumulator = SafeCalloc(dimensionCount, sizeof(double));
        for (size_t j = 0; j < stimulusCount; j++)
                cblas_daxpy((int)dimensionCount,
                            1.0 / (double)stimulusCount, 
                            normalisedCoordinates + dimensionCount * j, 
                            1, 
                            accumulator, 
                            1);
        for (size_t j = 0; j < stimulusCount; j++)
                cblas_daxpy((int)dimensionCount, 
                            -1.0, 
                            accumulator,
                            1, 
                            normalisedCoordinates + dimensionCount * j,
                            1);
        free(accumulator);
        if (weightMeans)
                for (size_t r = 0; r < dimensionCount; r++)
                        cblas_dscal((int)stimulusCount,
                                    sqrt(weightMeans[r]),
                                    normalisedCoordinates + r,
                                    (int)dimensionCount);
        return normalisedCoordinates;
}

/* N.B.: Basic initialisation must occur before calling. */
static double * 
NewNormalisedSpecificities(const ModelSpace * restrict self,
                           const double * restrict weightMeans,
                           const double * restrict specificities)
{
        if (!specificities) return NULL;
        double * restrict normalisedSpecificities;
        normalisedSpecificities = SafeMalloc(self->specificitiesSize,
                                             sizeof(double));
        cblas_dcopy((int)self->specificitiesSize,
                    specificities,
                    1,
                    normalisedSpecificities,
                    1);
        if (weightMeans
            && ModelSpecificityType(self->model) == GlobalSpecificities)
                cblas_dscal((int)StimulusCount(self->stimulusSet), 
                            weightMeans[DimensionCount(self->model)], 
                            normalisedSpecificities, 
                            1);
        return normalisedSpecificities;
}

static void InitialiseBasicValues(ModelSpace * self)
{
        if (self) {
                const size_t classCount = ClassCount(self->model);
                const size_t stimulusCount = StimulusCount(self->stimulusSet);
                const size_t pairCount = StimulusPairCount(self->stimulusSet);
                const size_t dimensionCount = DimensionCount(self->model);
                const SpecificityType sType = ModelSpecificityType(self->model);
                const size_t coordinatesSize = SizeProduct(stimulusCount, 
                                                           dimensionCount);
                if (!IsConvertibleToInt(coordinatesSize))
                        ExitWithError("Too many dimensions"
                                      " to compute a model space");
                size_t specificitiesSize;
                switch(sType) {
                        case ClassSpecificities:
                                specificitiesSize = SizeProduct(stimulusCount,
                                                                classCount);
                                break;
                        case GlobalSpecificities:
                                specificitiesSize = stimulusCount;
                                break;
                        default:
                                specificitiesSize = 0;
                                break;
                }
                const size_t differencesSize = SizeProduct(pairCount, 
                                                           dimensionCount);
                if (!IsConvertibleToInt(differencesSize))
                        ExitWithError("Too many dimensions"
                                      " to compute a model space");
                const size_t distancesSize = SizeProduct(pairCount, classCount);
                if (!IsConvertibleToInt(distancesSize))
                        ExitWithError("Too many classes"
                                      " to compute a model space");                
                self->coordinatesSize = coordinatesSize;
                self->specificitiesSize = specificitiesSize;
                self->differencesSize = differencesSize;
                self->distancesSize = distancesSize;
                self->parameterCount = ParameterCountValue(self);
        }
}

static void InitialiseDerivedValues(ModelSpace * self)
{
        if (self) {
                self->pairSpecificities = NewUnfoldedSpecificities(self);
                self->differences = NewDifferences(self);
                self->squaredDifferences = NewSquaredDifferences(self);
                self->distances = NewDistances(self);
        }
}

ModelSpace * NewModelSpace(const StimulusSet * restrict stimulusSet,
                           const Model * restrict model,
                           const double * restrict weights,
                           const double * restrict coordinates,
                           const double * restrict specificities)
{
        if (!stimulusSet || !model) return NULL;               
        ModelSpace * restrict self;
        if ((self = malloc(sizeof(ModelSpace)))) {
                self->stimulusSet = stimulusSet;
                self->model = model;
                InitialiseBasicValues(self);
                double * restrict wMeans = NULL;
                if (ClassCount(model) > 1) {
                        if (weights) {
                                wMeans = NewWeightMeans(self, weights);
                                self->weights = NewNormalisedWeights(self, 
                                                                     wMeans, 
                                                                     weights);                        
                        } else {
                                const size_t wSize = ExtendedWeightsSize(model);
                                self->weights = SafeMalloc(wSize, 
                                                           sizeof(double));
                                for (size_t w = 0; w < wSize; w++)
                                        self->weights[w] = 1.0;
                        }
                } else {
                        self->weights = NULL;
                }
                if (coordinates) {
                        self->coordinates 
                        = NewNormalisedCoordinates(self, 
                                                   wMeans, 
                                                   coordinates);
                } else {
                        double * restrict coords;
                        coords = SafeMalloc(self->coordinatesSize,
                                                 sizeof(double));
                        for (size_t x = 0; x < self->coordinatesSize; x++)
                                coords[x] = (-1.0 
                                             + (2.0 
                                                * (double)rand() 
                                                / (double)RAND_MAX));
                        self->coordinates = NewNormalisedCoordinates(self, 
                                                                     wMeans, 
                                                                     coords);
                }
                if (ModelSpecificityType(model)) {
                        if (specificities) {
                                self->specificities 
                                = NewNormalisedSpecificities(self,
                                                             wMeans,
                                                             specificities);
                        } else {
                                self->specificities 
                                = SafeMalloc(self->specificitiesSize, 
                                             sizeof(double));
                                for (size_t s = 0; 
                                     s < self->specificitiesSize; 
                                     s++)
                                        self->specificities[s] = 1.0;
                        }
                } else {
                        self->specificities = NULL;                        
                }
                free(wMeans);
                InitialiseDerivedValues(self);
        }
        return self;
}

ModelSpace * NewDummySpace(const StimulusSet * restrict stimulusSet,
                           const Model * restrict model,
                           double * restrict distances)
{
        if (!stimulusSet || !model || !distances) return NULL;               
        ModelSpace * restrict self;
        if ((self = malloc(sizeof(ModelSpace)))) {
                self->stimulusSet = stimulusSet;
                self->model = model;
                InitialiseBasicValues(self);
                self->weights = NULL;
                self->coordinates = NULL;
                self->specificities = NULL;                        
                self->pairSpecificities = NULL;
                self->differences = NULL;
                self->squaredDifferences = NULL;
                self->distances = distances;
        }
        return self;
}

ModelSpace * NewModelSpaceByUpdatingWeights(const ModelSpace * restrict self, 
                                            const double * restrict weights)
{
        if (!self || (!weights && self->weights)) return NULL;
        ModelSpace * restrict new;
        if ((new = malloc(sizeof(ModelSpace)))) {
                memcpy(new, self, sizeof(ModelSpace));
                if (ClassCount(self->model) > 1) {
                        const size_t wSize = ExtendedWeightsSize(new->model);
                        new->weights = SafeMalloc(wSize, sizeof(double));
                        cblas_dcopy((int)wSize, weights, 1, new->weights, 1);
                } else {
                        new->weights = NULL;
                }
                new->coordinates = SafeMalloc(self->coordinatesSize,
                                              sizeof(double));
                cblas_dcopy((int)self->coordinatesSize, 
                            self->coordinates, 
                            1, 
                            new->coordinates, 
                            1);
                if (ModelSpecificityType(self->model)) {
                        new->specificities = SafeMalloc((self
                                                          ->specificitiesSize),
                                                         sizeof(double));
                        cblas_dcopy((int)self->specificitiesSize, 
                                    self->specificities, 
                                    1, 
                                    new->specificities, 
                                    1);
                } else {
                        new->specificities = NULL;                        
                }
                InitialiseDerivedValues(new);
        }
        return new;
}

ModelSpace * NewModelSpaceByNormalisingWeights(const ModelSpace * restrict self)
{
        if (!self) return NULL;
        ModelSpace * restrict new;
        if ((new = malloc(sizeof(ModelSpace)))) {
                memcpy(new, self, sizeof(ModelSpace));
                double * restrict wMeans = NewWeightMeans(new, self->weights);
                new->weights = NewNormalisedWeights(new, wMeans, self->weights);
                new->coordinates = NewNormalisedCoordinates(new,
                                                            wMeans,
                                                            self->coordinates);
                new->specificities = NewNormalisedSpecificities(new,
                                                                wMeans,
                                                                self
                                                                ->specificities
                                                                );
                free(wMeans);
                InitialiseDerivedValues(new);
        }
        return new;
}

ModelSpace * 
NewModelSpaceByUpdatingCoordinates(const ModelSpace * restrict self, 
                                   const double * restrict coordinates,
                                   const double * restrict specificities)
{
        if (!self || !coordinates || (!specificities && self->specificities))
                return NULL;
        ModelSpace * restrict new;
        if ((new = malloc(sizeof(ModelSpace)))) {
                memcpy(new, self, sizeof(ModelSpace));
                if (new->weights) {
                        const size_t wSize = ExtendedWeightsSize(new->model);
                        new->weights = SafeMalloc(wSize, sizeof(double));
                        cblas_dcopy((int)wSize,
                                    self->weights,
                                    1,
                                    new->weights,
                                    1);                                                
                }
                new->coordinates = NewNormalisedCoordinates(new,
                                                            NULL,
                                                            coordinates);
                if (new->specificitiesSize) {
                        new->specificities = SafeMalloc(new->specificitiesSize,
                                                        sizeof(double));
                        cblas_dcopy((int)new->specificitiesSize, 
                                    specificities, 
                                    1, 
                                    new->specificities, 
                                    1);                        
                }
                InitialiseDerivedValues(new);
        }
        return new;
}

void DeleteModelSpace(ModelSpace * restrict self)
{
        if (self) {
                FreeAndClear(self->distances);
                FreeAndClear(self->squaredDifferences);
                FreeAndClear(self->differences);
                FreeAndClear(self->pairSpecificities);
                FreeAndClear(self->specificities);
                FreeAndClear(self->coordinates);
                FreeAndClear(self->weights);
                free(self);
        }
}

const StimulusSet * ModelSpaceStimulusSet(const ModelSpace * restrict self)
{
        return self ? self->stimulusSet : NULL;
}

const Model * ModelSpaceModel(const ModelSpace * restrict self)
{
        return self ? self->model : NULL;
}

const double * Weights(const ModelSpace * restrict self)
{
        return self ? self->weights : NULL;
}

const double * Coordinates(const ModelSpace * restrict self)
{
        return self ? self->coordinates : NULL;
}

const double * Specificities(const ModelSpace * restrict self)
{
        return self ? self->specificities : NULL;
}

size_t ParameterCount(const ModelSpace * restrict self)
{
        return self ? self->parameterCount : 0;
}

double DistanceForPairAndClass(const ModelSpace * restrict self,
                               const StimulusPair * pair,
                               size_t latentClass)
{
        if (self 
            && self->distances
            && pair 
            && latentClass < ClassCount(self->model)) {
                const size_t pairCount = StimulusPairCount(self->stimulusSet);
                const size_t number = PairNumberForStimulusPair(self
                                                                ->stimulusSet,
                                                                pair);
                return self->distances[pairCount * latentClass + number];
        }
        else return NAN;
}

const double * RawDifferences(const ModelSpace * self)
{
        return self ? self->differences : NULL;
}

const double * SquaredDifferences(const ModelSpace * self)
{
        return self ? self->squaredDifferences : NULL;
}

const double * PairwiseSpecificities(const ModelSpace * self)
{
        return self ? self->pairSpecificities : NULL;
}

const double * ClassDistances(const ModelSpace * self)
{
        return self ? self->distances : NULL;
}

size_t CoordinatesSize(const ModelSpace * self)
{
        return self ? self->coordinatesSize : 0;
}

size_t SpecificitiesSize(const ModelSpace * self)
{
        return self ? self->specificitiesSize : 0;
}

size_t DifferencesSize(const ModelSpace * self)
{
        return self ? self->differencesSize : 0;
}

size_t DistancesSize(const ModelSpace * self)
{
        return self ? self->distancesSize : 0;
}
