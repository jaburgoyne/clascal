// CLASCAL (CoordinateSolver.c)
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

#include <float.h>
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
#include "Experiment.h"
#include "Model.h"
#include "_Model.h"
#include "ClassAssignment.h"
#include "ModelSpace.h"
#include "_ModelSpace.h"
#include "Solution.h"
#include "_Solution.h"
#include "_CoordinateSolver.h"

typedef __CLPK_integer LPInteger;

struct _CoordinateSolver {
        Solution * solution; // non-const to allow for lazy initialisation
        size_t gradientSize;
        double * restrict gradient;
        size_t hessianSize;
        double * restrict hessian;
};
typedef struct _CoordinateSolver CoordinateSolver;

/*
 * Returns the gradient of q1 for the coordinates.
 */
static double * NewCoordinateGradient(const CoordinateSolver * restrict self)
{
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t coordinatesSize = CoordinatesSize(space);
        if (!coordinatesSize) return NULL;
        const double * restrict weightedRelativeErrors;
        weightedRelativeErrors = WeightedRelativeErrors(self->solution);
        if (!weightedRelativeErrors) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(stimulusSet);
        if (!stimulusPairs) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict weights = Weights(space);
        if (!weights && classCount > 1) return NULL;
        const double * restrict differences = RawDifferences(space);
        if (!differences) return NULL;
        const size_t weightsSize = WeightsSize(model);
        double * restrict gradient = SafeCalloc(coordinatesSize, 
                                                sizeof(double));
        double * restrict weightedDifferences = SafeMalloc(weightsSize, 
                                                           sizeof(double));
        double * restrict pairTerm = SafeCalloc(dimensionCount, sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                if (classCount > 1)
                        cblas_dcopy((int)weightsSize, 
                                    weights, 
                                    1, 
                                    weightedDifferences, 
                                    1);
                else 
                        for (size_t r = 0; r < weightsSize; r++)
                                weightedDifferences[r] = 1.0;
                for (size_t r = 0; r < dimensionCount; r++) {
                        cblas_dscal((int)classCount, 
                                    differences[dimensionCount * m + r], 
                                    weightedDifferences + classCount * r, 
                                    1);
                }
                cblas_dgemv(CblasRowMajor, 
                            CblasNoTrans, 
                            (int)dimensionCount, 
                            (int)classCount, 
                            -2.0, 
                            weightedDifferences, 
                            (int)classCount, 
                            weightedRelativeErrors + m, 
                            (int)pairCount, 
                            0.0, 
                            pairTerm, 
                            1);
                cblas_daxpy((int)dimensionCount, 
                            -1.0, 
                            pairTerm, 
                            1, 
                            gradient + dimensionCount * stimulusPairs[m].j, 
                            1);
                cblas_daxpy((int)dimensionCount, 
                            1.0, 
                            pairTerm, 
                            1, 
                            gradient + dimensionCount * stimulusPairs[m].k, 
                            1);
        }
        free(pairTerm);
        free(weightedDifferences);
        return gradient;
}

/*
 * Returns the gradient of q1 for the specificities.
 */
static double * NewSpecificityGradient(const CoordinateSolver * restrict self)
{
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t specificitiesSize = SpecificitiesSize(space);
        if (!specificitiesSize) return NULL;
        const double * restrict weightedRelativeErrors;
        weightedRelativeErrors = WeightedRelativeErrors(self->solution);
        if (!weightedRelativeErrors) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(space);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        if (!stimulusCount) return NULL;
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(stimulusSet);
        if (!stimulusPairs) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const SpecificityType specificityType = ModelSpecificityType(model);
        if (specificityType == NoSpecificities) return NULL;
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict weights = Weights(space);
        if (!weights && classCount > 1) return NULL;
        double * restrict gradient;
        gradient = SafeCalloc(specificitiesSize, sizeof(double));
        switch (specificityType) {
                case ClassSpecificities:
                        for (size_t m = 0; m < pairCount; m++) {
                                cblas_daxpy((int)classCount, 
                                            -1.0, 
                                            weightedRelativeErrors + m, 
                                            (int)pairCount, 
                                            gradient + stimulusPairs[m].j,
                                            (int)stimulusCount);                              
                                cblas_daxpy((int)classCount, 
                                            -1.0, 
                                            weightedRelativeErrors + m, 
                                            (int)pairCount, 
                                            gradient + stimulusPairs[m].k,
                                            (int)stimulusCount);                              
                        }
                        break;
                case GlobalSpecificities:
                        for (size_t m = 0; m < pairCount; m++) {
                                double dot = 0.0;
                                if (classCount > 1)
                                        dot = cblas_ddot((int)classCount, 
                                                         (weights 
                                                          + (classCount 
                                                             * dimensionCount)), 
                                                         1, 
                                                         (weightedRelativeErrors 
                                                          + m), 
                                                         (int)pairCount);
                                else dot = weightedRelativeErrors[m];
                                gradient[stimulusPairs[m].j] -= dot;
                                gradient[stimulusPairs[m].k] -= dot;
                        }
                        break;
                case NoSpecificities:
                        break;
        }
        return gradient;
}

/*
 * Returns the complete gradient of q1 with respect to dimensions and 
 * specificities. Appends the specificity gradient to the coordinate gradient.
 */
static double * NewGradient(const CoordinateSolver * restrict self)
{
        if (!self->gradientSize) return NULL;
        double * restrict coordinateGradient;
        coordinateGradient = NewCoordinateGradient(self);
        if (!coordinateGradient) return NULL;
        double * restrict specificityGradient;
        specificityGradient = NewSpecificityGradient(self);
        if (!specificityGradient) return coordinateGradient;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t coordinatesSize = CoordinatesSize(space);
        const size_t specificitiesSize = SpecificitiesSize(space);
        double * restrict gradient = SafeMalloc(self->gradientSize, 
                                                sizeof(double));
        cblas_dcopy((int)coordinatesSize,
                    coordinateGradient, 
                    1, 
                    gradient, 
                    1);
        cblas_dcopy((int)specificitiesSize,
                    specificityGradient, 
                    1, 
                    gradient + coordinatesSize, 
                    1);
        free(specificityGradient);
        free(coordinateGradient);
        return gradient;
}

/*
 * Returns the expected Hessian of coordinates against coordinates. The order
 * of the indices for rows and columns corresponds to that of the coordinate
 * gradient in row-major order. The matrix is unpacked. Note that Hessian
 * factors must be initialised before calling this function.
 */
static double * 
NewExpectedCoordinateHessian(const CoordinateSolver * restrict self)
{
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t coordinatesSize = CoordinatesSize(space);
        if (!coordinatesSize) return NULL;
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(stimulusSet);
        if (!stimulusPairs) return NULL;        
        const Model * restrict model = ModelSpaceModel(space);
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict weights = Weights(space);
        if (!weights && classCount > 1) return NULL;
        const double * restrict differences = RawDifferences(space);
        if (!differences) return NULL;
        const size_t hessianSize = SizeProduct(coordinatesSize, 
                                               coordinatesSize);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        double * restrict accumulator = SafeMalloc(classCount, sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                for (size_t r1 = 0; r1 < dimensionCount; r1++) {
                        for (size_t r2 = 0; r2 < dimensionCount; r2++) {
                                const double dr1 = differences[(dimensionCount 
                                                                * m) 
                                                               + r1];
                                const double dr2 = differences[(dimensionCount 
                                                                * m)
                                                               + r2];
                                if (classCount > 1) {
                                        cblas_dcopy((int)classCount, 
                                                    weights + classCount * r1, 
                                                    1, 
                                                    accumulator, 
                                                    1);
                                        cblas_dtbmv(CblasRowMajor, 
                                                    CblasUpper, 
                                                    CblasNoTrans, 
                                                    CblasNonUnit, 
                                                    (int)classCount, 
                                                    0, 
                                                    weights + classCount * r2, 
                                                    1, 
                                                    accumulator, 
                                                    1);
                                } else {
                                        for (size_t t = 0; t < classCount; t++)
                                                accumulator[t] = 1.0;
                                }
                                const double dot = cblas_ddot((int)classCount, 
                                                              (hessianFactors 
                                                               + m), 
                                                              (int)pairCount, 
                                                              accumulator, 
                                                              1);
                                const double a = 2.0 * dr1 * dr2 * dot;
                                const size_t jr1 = ((dimensionCount 
                                                     * stimulusPairs[m].j)
                                                    + r1);
                                const size_t jr2 = ((dimensionCount 
                                                     * stimulusPairs[m].j)
                                                    + r2);
                                const size_t kr1 = ((dimensionCount 
                                                     * stimulusPairs[m].k)
                                                    + r1);
                                const size_t kr2 = ((dimensionCount 
                                                     * stimulusPairs[m].k)
                                                    + r2);
                                hessian[coordinatesSize * kr1 + jr2] -= a;
                                hessian[coordinatesSize * jr2 + kr1] -= a;
                                hessian[coordinatesSize * jr1 + jr2] += a;
                                hessian[coordinatesSize * kr1 + kr2] += a;
                        }
                }
        }
        free(accumulator);
        return hessian;
}

/*
 * Returns the expected Hessian of coordinates against specificities. The order
 * of the indices for rows corresponds to that of the coordinate gradient in 
 * row-major order; for the columns, to that of the specificity gradient. The 
 * matrix is unpacked.  Note that Hessian factors must be initialised before 
 * calling this function.
 */
static double * 
NewExpectedCoordinateSpecificityHessian(const CoordinateSolver * restrict self)
{
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t coordinatesSize = CoordinatesSize(space);
        if (!coordinatesSize) return NULL;
        const size_t specificitiesSize = SpecificitiesSize(space);
        if (!specificitiesSize) return NULL;
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(stimulusSet);
        if (!stimulusPairs) return NULL;        
        const Model * restrict model = ModelSpaceModel(space);
        const SpecificityType specificityType = ModelSpecificityType(model);
        if (specificityType == NoSpecificities) return NULL;
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict weights = Weights(space);
        if (!weights && classCount > 1) return NULL;
        const double * restrict differences = RawDifferences(space);
        if (!differences) return NULL;
        const size_t hessianSize = SizeProduct(coordinatesSize, 
                                               specificitiesSize);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                for (size_t t = 0; t < classCount; t++) {
                        for (size_t r = 0; r < dimensionCount; r++) {
                                const double h = hessianFactors[pairCount * t 
                                                                + m];
                                const double wr = (classCount > 1
                                                   ? weights[classCount * r + t]
                                                   : 1.0);
                                const double ws = (((classCount > 1)
                                                    && (specificityType 
                                                        == GlobalSpecificities))
                                                   ? weights[(classCount 
                                                              * dimensionCount)
                                                             + t]
                                                   : 1.0);
                                const double d = differences[dimensionCount * m
                                                             + r];
                                const double term = h * wr * ws * d;
                                const size_t xjr = ((dimensionCount 
                                                     * stimulusPairs[m].j) 
                                                    + r);
                                const size_t xkr = ((dimensionCount 
                                                     * stimulusPairs[m].k)
                                                    + r);
                                const size_t stj = ((specificityType 
                                                     == ClassSpecificities)
                                                    ? (stimulusCount * t 
                                                       + stimulusPairs[m].j)
                                                    : stimulusPairs[m].j);
                                const size_t stk = ((specificityType 
                                                     == ClassSpecificities)
                                                    ? (stimulusCount * t
                                                       + stimulusPairs[m].k)
                                                    : stimulusPairs[m].k);
                                hessian[specificitiesSize * xjr + stj] -= term;
                                hessian[specificitiesSize * xkr + stj] += term;
                                hessian[specificitiesSize * xjr + stk] -= term;
                                hessian[specificitiesSize * xkr + stk] += term;
                        }
                }
        }
        return hessian;
}

/*
 * Returns the expected Hessian of specificities against specificities. The 
 * order of the indices corresponds to that of the specificity gradient. The 
 * matrix is unpacked.  Note that Hessian factors must be initialised before 
 * calling this function.
 */
static double * 
NewExpectedSpecificityHessian(const CoordinateSolver * restrict self)
{
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t specificitiesSize = SpecificitiesSize(space);
        if (!specificitiesSize) return NULL;
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(space);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const StimulusPair * restrict stimulusPairs;
        stimulusPairs = StimulusPairs(stimulusSet);
        if (!stimulusPairs) return NULL;        
        const Model * restrict model = ModelSpaceModel(space);
        const SpecificityType specificityType = ModelSpecificityType(model);
        if (specificityType == NoSpecificities) return NULL;
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict weights = Weights(space);
        if (!weights && classCount > 1) return NULL;
        const size_t hessianSize = SizeProduct(specificitiesSize, 
                                               specificitiesSize);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        for (size_t m = 0; m < pairCount; m++) {
                for (size_t t = 0; t < classCount; t++) {
                        const double h = hessianFactors[pairCount * t + m];
                        const double w = (((classCount > 1)
                                           && (specificityType 
                                               == GlobalSpecificities))
                                          ? weights[classCount * dimensionCount 
                                                    + t]
                                          : 1.0);
                        const double term = 0.5 * h * w * w;
                        const size_t stj = ((specificityType 
                                             == ClassSpecificities)
                                            ? (stimulusCount * t + 
                                               stimulusPairs[m].j)
                                            : stimulusPairs[m].j);
                        const size_t stk = ((specificityType 
                                             == ClassSpecificities)
                                            ? (stimulusCount * t 
                                               + stimulusPairs[m].k)
                                            : stimulusPairs[m].k);
                        hessian[specificitiesSize * stj + stj] += term;
                        hessian[specificitiesSize * stj + stk] += term;
                        hessian[specificitiesSize * stk + stj] += term;
                        hessian[specificitiesSize * stk + stk] += term;
                }
        }
        return hessian;        
}

static CoordinateSolver * 
NewCoordinateSolver(Solution * restrict solution)
{
        if (!solution) return NULL;
        const ModelSpace * space = SolutionModelSpace(solution);
        size_t coordinatesSize = CoordinatesSize(space);
        if (!coordinatesSize) return NULL;
        size_t specificitiesSize = SpecificitiesSize(space);
        size_t gradientSize = SizeSum(coordinatesSize, specificitiesSize);
        size_t hessianSize = SizeProduct(gradientSize, gradientSize);
        if (!IsConvertibleToInt(hessianSize))
                ExitWithError("Model too large to solve for coordinates");
        CoordinateSolver * restrict self;
        if ((self = malloc(sizeof(CoordinateSolver)))) {
                self->solution = solution;
                self->gradientSize = gradientSize;
                self->gradient = NewGradient(self);
                self->hessianSize = hessianSize;
                self->hessian = NULL;
        }
        return self;
}

static void DeleteCoordinateSolver(CoordinateSolver * restrict self)
{
        if (self) {
                FreeAndClear(self->hessian);
                FreeAndClear(self->gradient);
                free(self);
        }
}

/*
 * Returns the full expected Hessian, with indices corresponding to the
 * full gradient. The matrix is unpacked.
 */
static double * ExpectedHessian(CoordinateSolver * restrict self)
{
        if (self->hessian) return self->hessian;
        double * restrict coordinateHessian;
        coordinateHessian = NewExpectedCoordinateHessian(self);
        if (!coordinateHessian) return NULL;
        double * restrict jointHessian;
        jointHessian = NewExpectedCoordinateSpecificityHessian(self);
        if (!jointHessian) return coordinateHessian;
        double * restrict specificityHessian;
        specificityHessian = NewExpectedSpecificityHessian(self);
        if (!specificityHessian) {
                free(jointHessian);
                return coordinateHessian;
        }
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const size_t coordinatesSize = CoordinatesSize(space);
        const size_t specificitiesSize = SpecificitiesSize(space);
        self->hessian = SafeMalloc(self->hessianSize, sizeof(double));
        for (size_t x = 0; x < coordinatesSize; x++) {
                cblas_dcopy((int)coordinatesSize, 
                            coordinateHessian + coordinatesSize * x, 
                            1, 
                            self->hessian + self->gradientSize * x, 
                            1);
                cblas_dcopy((int)specificitiesSize,
                            jointHessian + specificitiesSize * x,
                            1,
                            (self->hessian + 
                             self->gradientSize * x 
                             + coordinatesSize),
                            1);
        }
        for (size_t s = 0; s < specificitiesSize; s++) {
                cblas_dcopy((int)coordinatesSize,
                            jointHessian + s,
                            (int)specificitiesSize,
                            (self->hessian 
                             + self->gradientSize * (coordinatesSize + s)),
                            1);
                cblas_dcopy((int)specificitiesSize,
                            specificityHessian + specificitiesSize * s,
                            1,
                            (self->hessian 
                             + self->gradientSize * (coordinatesSize + s)
                             + coordinatesSize),
                            1);
        }
        free(specificityHessian);
        free(jointHessian);
        free(coordinateHessian);
        return self->hessian;
}

/**,
 * Returns true if the gradient is negative for a frozen specificity.
 */
static bool FrozenGradientIsNegative(const CoordinateSolver * restrict self,
                                     Solution * restrict solution0)
{
        if (!self || !self->gradient) return false;
        const Parameters * restrict parameters;
        parameters = SolutionParameters(self->solution);
        const ModelSpace * restrict space = SolutionModelSpace(solution0);
        const double * restrict specificities = Specificities(space);
#ifdef WINSBERG_LEGACY
        // Winsberg's code handles class specificities incorrectly.
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        const size_t dimensionCount = DimensionCount(ModelSpaceModel(space));
        for (size_t s = 0; s < stimulusCount; s++)
#else
        const size_t coordinatesSize = CoordinatesSize(space);
        const size_t specificitiesSize = SpecificitiesSize(space);
        for (size_t s = 0; s < specificitiesSize; s++)
#endif
                if (isless(specificities[s], 
                           parameters->minSpecificityValue + parameters->margin)
#ifdef WINSBERG_LEGACY
                    // Winsberg's code misaddresses the gradient.
                    && isless(self->gradient[dimensionCount * s], -1e-3))
#else
                    && isless(self->gradient[coordinatesSize + s], 0.0))
#endif
                        return true;
        return false;
}

/**
 * Returns the search direction for a line search for coordinates: in this case,
 * the (Moore-Penrose pseudo)-inverse Hessian times the gradient.
 */
static double * NewSearchDirection(CoordinateSolver * restrict self,
                                   Solution * restrict solution0)
{
        if (!self || !self->gradient) return NULL;
        double * hessian = ExpectedHessian(self);
        if (!hessian) return NULL;
        const Parameters * restrict params = SolutionParameters(solution0);
        const ModelSpace * restrict space = SolutionModelSpace(solution0);
        const size_t coordinatesSize = CoordinatesSize(space);
        const double * restrict specificities = Specificities(space);
        size_t reducedGradientSize = self->gradientSize;
        // The following allocation is possibly too large, but it allows the for
        // loop following to fill the permutation at the same time as it 
        // determines the reduced size.
        size_t * restrict permutation = SafeCalloc(self->gradientSize,
                                                   sizeof(size_t));
        for (size_t x = 0; x < coordinatesSize; x++) permutation[x] = x;
        // The following loop will not run if there are no specificities.
        size_t nextIndex = coordinatesSize;
        double minGradValue = 0.0;
        size_t minGradIndex = self->gradientSize;
        for (size_t s = 0; s < SpecificitiesSize(space); s++) {
                if (isless(specificities[s], 
                           (params->minSpecificityValue + params->margin))) {
                        reducedGradientSize--;
                        if (isless(self->gradient[coordinatesSize + s],
                                   minGradValue)) {
                                minGradIndex = coordinatesSize + s;
                                minGradValue = self->gradient[minGradIndex];
                        }
                }
                else permutation[nextIndex++] = coordinatesSize + s;                
        }
#ifndef WINSBERG_LEGACY
        // Winsberg has a bug that effectively prevents specificities from ever
        // being unfrozen.
        if (isless(minGradValue, 0.0) 
            // Check the gradient size just in case something went wrong.
            && reducedGradientSize < self->gradientSize) {
                permutation[nextIndex] = minGradIndex;
                reducedGradientSize++;
        }
#endif
        const size_t reducedHessianSize = SizeProduct(reducedGradientSize,
                                                      reducedGradientSize);
        double * restrict invReducedGradient;
        invReducedGradient = SafeMalloc(reducedGradientSize, sizeof(double));
        double * restrict reducedHessian;
        reducedHessian = SafeMalloc(reducedHessianSize, sizeof(double));
        if (reducedGradientSize == self->gradientSize) {
                cblas_dcopy((int)reducedGradientSize, 
                            self->gradient, 
                            1, 
                            invReducedGradient, 
                            1);
                cblas_dscal((int)reducedGradientSize, 
                            -1.0, 
                            invReducedGradient, 
                            1);
                cblas_dcopy((int)reducedHessianSize,
                            hessian,
                            1,
                            reducedHessian,
                            1);
        } else {
                for (size_t a = 0; a < reducedGradientSize; a++) {
                        invReducedGradient[a] = -self->gradient[permutation[a]];
                        for (size_t b = 0; b < reducedGradientSize; b++)
                                reducedHessian[reducedGradientSize * a + b] 
                                = hessian[self->gradientSize * permutation[a] 
                                          + permutation[b]];
                }
        }
        // LAPACK is all call-by-reference.
        LPInteger m = (LPInteger)reducedGradientSize;
        LPInteger n = (LPInteger)reducedGradientSize;
        LPInteger nrhs = 1;
        LPInteger lda = (LPInteger)reducedGradientSize;
        LPInteger ldb = (LPInteger)reducedGradientSize;
        double * s = SafeMalloc(reducedGradientSize, sizeof(double));
        double rcond = params->invConditionNumber;
        LPInteger rank;
        double * work = SafeMalloc(1, sizeof(double));
        LPInteger lwork = -1;
        LPInteger info = 0;
        // This call determines work size only.
        // DGELSD refuses to compute LIWORK on OS X, so we use DGELSS. 
        dgelss_(&m, 
                &n, 
                &nrhs, 
                reducedHessian, 
                &lda, 
                invReducedGradient, 
                &ldb, 
                s, 
                &rcond, 
                &rank, 
                work, 
                &lwork, 
                &info);
        lwork = (LPInteger)lround(work[0]);
        free(work);
        work = SafeMalloc((size_t)lwork, sizeof(double));
        // This call is the real thing.
        dgelss_(&m, 
                &n, 
                &nrhs, 
                reducedHessian, 
                &lda, 
                invReducedGradient, 
                &ldb, 
                s, 
                &rcond, 
                &rank, 
                work, 
                &lwork, 
                &info);
        free(work); 
        free(s);
        double * restrict searchDirection;
        searchDirection = SafeCalloc(self->gradientSize, sizeof(double));
        for (size_t a = 0; a < reducedGradientSize; a++)
                // LAPACK has overwritten invReducedGradient with the solution.
                searchDirection[permutation[a]] = invReducedGradient[a];
        free(reducedHessian);
        free(invReducedGradient);
        free(permutation);
        return searchDirection;
}

// This routine performs the line search. It follows Winsberg directly.
// N.B.: It may return the same solution as the one passed in.
// WARNING: This is a simple line search that relies on starting with a
//          negative slope. With an insufficient number of iterations, it can
//          return a solution with positive slope.
#ifdef WINSBERG_LEGACY
static Solution * NewLineSearchSolution(CoordinateSolver * restrict self,
                                        Solution * restrict solution0,
                                        double SSR)
#else
static Solution * NewLineSearchSolution(CoordinateSolver * restrict self,
                                        Solution * restrict solution0)
#endif
{
        if (!solution0) return NULL;
        const Parameters * restrict parameters = SolutionParameters(solution0);
        double * restrict searchDirection0;
        searchDirection0 = NewSearchDirection(self, solution0);
        if (!searchDirection0) return NULL;
        const double norm = cblas_dnrm2((int)self->gradientSize, 
                                        searchDirection0, 
                                        1);
        cblas_dscal((int)self->gradientSize,
                    1.0 / norm,
                    searchDirection0,
                    1);
        const double slope = cblas_ddot((int)self->gradientSize,
                                        self->gradient,
                                        1,
                                        searchDirection0,
                                        1);
        if (isgreaterequal(slope, 0.0)) 
                // ExitWithError("Optimisation slope is non-negative");
                return solution0;
#ifndef WINSBERG_LEGACY
        const double SSR = SumOfSquaredModelError(self->solution);
#endif
        if (parameters->verbosity >= VERY_VERY_VERY_VERBOSE)
                fprintf(stdout, 
                        "                        \n"
                        "                        Line search\n"
                        "                        \n"
                        "                        Squared model error after"
                        " iteration 0: %e\n",
                        SSR);
        double stepSizes[5] = {0.0, 0.0, 0.0, NAN, NAN};
        double errors[5] = {SSR, SSR, SSR, NAN, NAN};
        double slopes[5] = {slope, slope, slope, NAN, NAN};
        stepSizes[4] = parameters->firstStepSize;
        const ModelSpace * restrict space0 = SolutionModelSpace(solution0);
        const double * restrict coordinates = Coordinates(space0);
        const size_t coordinatesSize = CoordinatesSize(space0);
        const double * restrict specificities = Specificities(space0);
        const size_t specificitiesSize = SpecificitiesSize(space0);
        const Model * restrict model = ModelSpaceModel(space0);
        const SpecificityType specificityType = ModelSpecificityType(model);
        double maxStepSize = INFINITY;
        if (specificityType) {
                for (size_t s = 0; s < specificitiesSize; s++) {
                        const size_t i = coordinatesSize + s;
                        if (isgreaterequal(searchDirection0[i], 0.0)) continue;
                        const double thisMax = ((parameters->minSpecificityValue
                                                 - specificities[s])
                                                / searchDirection0[i]);
                        maxStepSize = fmin(maxStepSize, thisMax);
                }
        }
        if (islessequal(maxStepSize, 0.0)) return solution0;
        stepSizes[4] = fmin(stepSizes[4], maxStepSize);
        Solution * restrict solution = solution0;
        CoordinateSolver * restrict solver = NULL;
        double * restrict nextCoordinates = SafeMalloc(coordinatesSize,
                                                       sizeof(double));
        double * restrict nextSpecs = (specificityType
                                       ? SafeMalloc(specificitiesSize,
                                                    sizeof(double))
                                       : NULL);
        for (size_t i = 0; i < parameters->maxSearchIterationCount; i++) {
                cblas_dcopy((int)coordinatesSize, 
                            coordinates, 
                            1, 
                            nextCoordinates, 
                            1);
                cblas_daxpy((int)coordinatesSize, 
                            stepSizes[4], 
                            searchDirection0, 
                            1, 
                            nextCoordinates,
                            1);
                if (specificityType) {
                        cblas_dcopy((int)specificitiesSize, 
                                    specificities, 
                                    1, 
                                    nextSpecs, 
                                    1);
                        cblas_daxpy((int)specificitiesSize, 
                                    stepSizes[4], 
                                    searchDirection0 + coordinatesSize, 
                                    1, 
                                    nextSpecs,
                                    1);
                }
                Solution * restrict newSolution;
                newSolution = NewSolutionByUpdatingCoordinates(solution0, 
                                                               nextCoordinates, 
                                                               nextSpecs);
                if (solution != solution0)
                        DeleteSolutionPreservingClassAssignment(solution);
                solution = newSolution;
                errors[4] = SumOfSquaredModelError(solution);
                if (parameters->verbosity >= VERY_VERY_VERY_VERBOSE)
                        fprintf(stdout, 
                                "                        Squared model error"
                                " after iteration %zu: %e\n",
                                SizeSum(i, 1),
                                errors[4]);
                DeleteCoordinateSolver(solver);
                solver = NewCoordinateSolver(solution);
                slopes[4] = cblas_ddot((int)solver->gradientSize, 
                                       solver->gradient,
                                       1,
                                       searchDirection0,
                                       1);
                const bool hasWolfeCurvature = isless(fabs(slopes[4]),
                                                      (parameters
                                                       ->slopeImprovementFactor)
                                                      * fabs(slopes[0]));
                const bool errorIsWorse = isgreater(errors[4], errors[0]);
                const bool slopeIsPositive = isgreater(slopes[4], 0.0);
                if (errorIsWorse && (hasWolfeCurvature || !slopeIsPositive)) {
                        // Backtracking
                        stepSizes[4] *= parameters->stepDecreaseFactor;
                        stepSizes[2] = stepSizes[1] = stepSizes[0];
                        errors[2] = errors[1] = errors[0];
                        slopes[2] = slopes[1] = slopes[0];
                } else {
                        if (hasWolfeCurvature) break;
                        if (!slopeIsPositive) { // Non-neg slope & better error
                                stepSizes[4] *= parameters->stepIncreaseFactor;
                        } else { // Positive slope regardless of error.
                                stepSizes[3] = stepSizes[4];
                                errors[3] = errors[4];
                                slopes[3] = slopes[4];
                                double z = (((3.0 
                                              / (stepSizes[4] - stepSizes[2]))
                                             * (errors[2] - errors[4]))
                                            + slopes[2]
                                            + slopes[4]);
                                double w1 = z * z - slopes[2] * slopes[4];
                                if (isless(fabs(slopes[2] 
                                                + slopes[4] 
                                                + 2.0 * z),
                                           1e-5)
                                    || isless(w1, 0.0)) {
                                        stepSizes[4] = (stepSizes[2]
                                                        - (slopes[2]
                                                           * ((stepSizes[4]
                                                               - stepSizes[2])
                                                              / (slopes[4]
                                                                 - slopes[2])))
                                                        );
                                } else {
                                        w1 = sqrt(w1);
                                        stepSizes[4] = (stepSizes[2]
                                                        + ((1.0
                                                            - ((slopes[4]
                                                                + w1
                                                                - z)
                                                               / (slopes[4]
                                                                  - slopes[2]
                                                                  + 2.0 * w1)))
                                                           * (stepSizes[4]
                                                              - stepSizes[2])));
                                }
                        }
                }
                if (islessequal(stepSizes[4], 0.0)) break; // convergence
                if (specificityType) {
                        if (stepSizes[4] 
                            == maxStepSize * parameters->stepIncreaseFactor)
                                break;
                        stepSizes[4] = fmin(stepSizes[4], maxStepSize);
                }
        }
        if (parameters->verbosity >= VERY_VERY_VERY_VERBOSE) 
                fprintf(stdout, "\n");
        free(nextCoordinates);
        free(nextSpecs);
        DeleteCoordinateSolver(solver);
        free(searchDirection0);
        return solution;
}

Solution * NewSpatialSolution(Solution * restrict solution0)
{
        if (!solution0) return NULL;
        const Parameters * restrict params = SolutionParameters(solution0);
        CoordinateSolver * restrict solver = NewCoordinateSolver(solution0);
        Solution * restrict soln = solution0;
        double SSR = SumOfSquaredModelError(soln);
        if (params->verbosity >= VERY_VERY_VERBOSE)
                fprintf(stdout,
                        "                \n"
                        "                Optimisation of coordinates\n"
                        "                \n"
                        "                Squared model error after iteration"
                        " 0: %e\n",
                        SSR);
#ifdef WINSBERG_LEGACY
        CoordinateSolver * restrict solver0 = solver;
        size_t iccount = 0;
#endif
        for (size_t i = 0; i < params->maxSpatialIterationCount; i++) {
#ifdef WINSBERG_LEGACY
        Winsberg:
                if (++iccount > 7) break;
#endif
                Solution * restrict newSolution;
#ifdef WINSBERG_LEGACY
                // For the first iteration only (and repeats of it due to goto),
                // the solver is not overwritten by solver0 in Winsberg's code.
                if (i) {
                        newSolution = NewLineSearchSolution(solver0, soln, SSR);                        
                } else {
                        newSolution = NewLineSearchSolution(solver, soln, SSR);
                }
#else
                newSolution = NewLineSearchSolution(solver, soln);
#endif
                if (newSolution == soln) break; // convergence
                double newSSR = SumOfSquaredModelError(newSolution);
                if (params->verbosity >= VERY_VERY_VERBOSE)
                        fprintf(stdout,
                                "                Squared model error after"
                                " iteration %zu: %e\n",
                                SizeSum(i, 1),
                                newSSR);
                double relativeImprovement = (SSR - newSSR) / SSR;
                if (isgreater(relativeImprovement, 0.0)) {
                        if (soln != solution0) 
                                DeleteSolutionPreservingClassAssignment(soln);
                        soln = newSolution;
#ifdef WINSBERG_LEGACY
                        if (solver != solver0) DeleteCoordinateSolver(solver);
#else
                        DeleteCoordinateSolver(solver);
#endif
                        solver = NewCoordinateSolver(soln);
                        if (SpecificitiesSize(SolutionModelSpace(soln))
                            && isless(relativeImprovement, 
                                      params->spatialImprovementFactor)) {
                                bool isAtBound;
                                isAtBound = FrozenGradientIsNegative(solver,
                                                                     soln);
#ifdef WINSBERG_LEGACY
                                if (isAtBound) goto Winsberg; else break;
#else
                                if (!isAtBound) break; 
#endif
                        }
                } else {
                        // We would have hit break above if newSln == solution0
                        DeleteSolutionPreservingClassAssignment(newSolution);
                        break;
                }
#ifndef WINSBERG_LEGACY
                SSR = newSSR;
#endif
        }
        if (params->verbosity >= VERY_VERY_VERBOSE) fprintf(stdout, "\n");
#ifdef WINSBERG_LEGACY
        if (solver0 != solver) DeleteCoordinateSolver(solver0);
#endif
        DeleteCoordinateSolver(solver);
        return soln;
}
