// CLASCAL (WeightSolver.c)
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
#include "_WeightSolver.h"

typedef __CLPK_integer LPInteger;

struct _WeightSolver {
        Solution * solution; // non-const to allow for lazy initialisation
        size_t gradientSize;
        double * restrict gradient;
        size_t hessianSize;
        double * restrict hessian;
};
typedef struct _WeightSolver WeightSolver;

/*
 * Returns the gradient of q1 for the weights.
 */
static double * NewGradient(const WeightSolver * restrict self)
{
        const double * restrict weightedRelativeErrors;
        weightedRelativeErrors = WeightedRelativeErrors(self->solution);
        if (!weightedRelativeErrors) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const double * restrict squaredDifferences = SquaredDifferences(space);
        if (!squaredDifferences) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL;
        const SpecificityType specificityType = ModelSpecificityType(model);
        const double * restrict pairwiseSpecificities;
        pairwiseSpecificities = PairwiseSpecificities(space);
        if (specificityType == GlobalSpecificities && !pairwiseSpecificities) 
                return NULL;
        const size_t weightsSize = WeightsSize(model);
        if (specificityType == GlobalSpecificities && !weightsSize) return NULL;
        double * restrict gradient = SafeCalloc(self->gradientSize, 
                                                sizeof(double));
        cblas_dgemm(CblasRowMajor, 
                    CblasTrans, 
                    CblasTrans, 
                    (int)dimensionCount, 
                    (int)classCount, 
                    (int)pairCount, 
                    -1.0, 
                    squaredDifferences, 
                    (int)dimensionCount, 
                    weightedRelativeErrors, 
                    (int)pairCount, 
                    0.0, 
                    gradient,
                    (int)classCount);
        if (specificityType == GlobalSpecificities)
                cblas_dgemv(CblasRowMajor, 
                            CblasNoTrans, 
                            (int)classCount,
                            (int)pairCount, 
                            -1.0, 
                            weightedRelativeErrors, 
                            (int)pairCount, 
                            pairwiseSpecificities, 
                            1, 
                            0.0, 
                            gradient + weightsSize, 
                            1);
        return gradient;
}

/*
 * Returns the expected Hessian of coordinate weights against coordinate 
 * weights. The order of the indices for rows and columns corresponds to that 
 * of the coordinate gradient in row-major order. The matrix is unpacked. Note 
 * that Hessian factors must be initialised before calling this function.
 */
static double * 
NewExpectedCoordinateHessian(const WeightSolver * restrict self)
{
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const double * restrict squaredDifferences = SquaredDifferences(space);
        if (!squaredDifferences) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL;
        const size_t weightsSize = WeightsSize(model);
        if (!weightsSize) return NULL;
        const size_t hessianSize = SizeProduct(weightsSize, weightsSize);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        double * restrict accumulator = SafeMalloc(pairCount, sizeof(double));
        for (size_t t = 0; t < classCount; t++) {
                for (size_t r1 = 0; r1 < dimensionCount; r1++) {
                        for (size_t r2 = 0; r2 < dimensionCount; r2++) {                                
                                cblas_dcopy((int)pairCount,
                                            squaredDifferences + r1,
                                            (int)dimensionCount,
                                            accumulator,
                                            1);
                                cblas_dtbmv(CblasRowMajor, 
                                            CblasUpper, 
                                            CblasNoTrans, 
                                            CblasNonUnit, 
                                            (int)pairCount, 
                                            0,
                                            squaredDifferences + r2, 
                                            (int)dimensionCount, 
                                            accumulator,
                                            1); // elementwise multiplication
                                const double h = (0.5 
                                                  * cblas_ddot((int)pairCount,
                                                               (hessianFactors
                                                                + (pairCount 
                                                                   * t)),
                                                               1,
                                                               accumulator,
                                                               1));
                                const size_t tr1 = classCount * r1 + t;
                                const size_t tr2 = classCount * r2 + t;
#ifdef WINSBERG_LEGACY
                                // Winsberg inadvertently doubles the values of
                                // coordinate-coordinate Hessian entries because
                                // she effectively uses += rather than =.
                                hessian[weightsSize * tr1 + tr2] += h;
                                hessian[weightsSize * tr2 + tr1] += h;                                
#else
                                hessian[weightsSize * tr1 + tr2] = h;
                                hessian[weightsSize * tr2 + tr1] = h;
#endif
                        }
                }
        }
        FreeAndClear(accumulator);
        return hessian;
}

/*
 * Returns the expected Hessian of coordinate weights against specificity 
 * weights. The order of the indices for rows and columns corresponds to that 
 * of the coordinate gradient in row-major order. The matrix is unpacked. Note 
 * that Hessian factors must be initialised before calling this function.
 */
static double * 
NewExpectedCoordinateSpecificityHessian(const WeightSolver * restrict self)
{
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const double * restrict squaredDifferences = SquaredDifferences(space);
        if (!squaredDifferences) return NULL;
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const size_t dimensionCount = DimensionCount(model);
        if (!dimensionCount) return NULL;
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL;
        const size_t weightsSize = WeightsSize(model);
        if (!weightsSize) return NULL;
        const SpecificityType specificityType = ModelSpecificityType(model);
        if (specificityType != GlobalSpecificities) return NULL;
        const double * restrict pairwiseSpecificities;
        pairwiseSpecificities = PairwiseSpecificities(space);
        if (!pairwiseSpecificities) return NULL;
        const size_t hessianSize = SizeProduct(weightsSize, classCount);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        double * restrict accumulator = SafeMalloc(pairCount, sizeof(double));
        for (size_t t = 0; t < classCount; t++) {
                for (size_t r = 0; r < dimensionCount; r++) {
                        cblas_dcopy((int)pairCount,
                                    squaredDifferences + r,
                                    (int)dimensionCount,
                                    accumulator,
                                    1);
                        cblas_dtbmv(CblasRowMajor, 
                                    CblasUpper, 
                                    CblasNoTrans, 
                                    CblasNonUnit, 
                                    (int)pairCount, 
                                    0,
                                    pairwiseSpecificities, 
                                    1, 
                                    accumulator,
                                    1); // trick for elementwise multiplication
                        const double h = 0.5 * cblas_ddot((int)pairCount,
                                                          (hessianFactors
                                                           + pairCount * t),
                                                          1,
                                                          accumulator,
                                                          1);
                        const size_t tr = classCount * r + t;
                        hessian[classCount * tr + t] = h;
                }
        }
        FreeAndClear(accumulator);
        return hessian;
}

/*
 * Returns the expected Hessian of coordinate weights against specificity 
 * weights. The order of the indices for rows and columns corresponds to that 
 * of the coordinate gradient in row-major order. The matrix is unpacked. Note 
 * that Hessian factors must be initialised before calling this function.
 */
static double * 
NewExpectedSpecificityHessian(const WeightSolver * restrict self)
{
        const double * restrict hessianFactors = HessianFactors(self->solution);
        if (!hessianFactors) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const StimulusSet * restrict stimulusSet = ModelSpaceStimulusSet(space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const Model * restrict model = ModelSpaceModel(space);
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL;
        const SpecificityType specificityType = ModelSpecificityType(model);
        if (specificityType != GlobalSpecificities) return NULL;
        const double * restrict pairwiseSpecificities;
        pairwiseSpecificities = PairwiseSpecificities(space);
        if (!pairwiseSpecificities) return NULL;
        const size_t hessianSize = SizeProduct(classCount, classCount);
        double * restrict hessian = SafeCalloc(hessianSize, sizeof(double));
        double * restrict squaredPairwiseSpecificities;
        squaredPairwiseSpecificities = SafeCalloc(pairCount, sizeof(double));
        double * restrict accumulator = SafeCalloc(classCount, sizeof(double));
        cblas_dcopy((int)pairCount,
                    pairwiseSpecificities,
                    1,
                    squaredPairwiseSpecificities,
                    1);
        cblas_dtbmv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)pairCount, 
                    0,
                    pairwiseSpecificities, 
                    1, 
                    squaredPairwiseSpecificities,
                    1); // trick for elementwise multiplication
        cblas_dgemv(CblasRowMajor, 
                    CblasNoTrans, 
                    (int)classCount,
                    (int)pairCount, 
                    0.5, 
                    hessianFactors, 
                    (int)pairCount, 
                    squaredPairwiseSpecificities, 
                    1, 
                    0.0, 
                    accumulator, 
                    1);
        for (size_t t = 0; t < classCount; t++) 
                hessian[classCount * t + t] = accumulator[t];
        FreeAndClear(accumulator);
        FreeAndClear(squaredPairwiseSpecificities);
        return hessian;
}

static WeightSolver * NewWeightSolver(Solution * restrict solution)
{
        if (!solution) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(solution);
        const Model * restrict model = ModelSpaceModel(space);
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL; // weights are all 1 for a single class
        const size_t gradientSize = ExtendedWeightsSize(model);
        size_t hessianSize = SizeProduct(gradientSize, gradientSize);
        if (!IsConvertibleToInt(hessianSize))
                ExitWithError("Model too large to solve for weights");
        WeightSolver * restrict self;
        if ((self = malloc(sizeof(WeightSolver)))) {
                self->solution = solution;
                self->gradientSize = gradientSize;
                self->gradient = NewGradient(self);
                self->hessianSize = hessianSize;
                self->hessian = NULL;
        }
        return self;
}

static void DeleteWeightSolver(WeightSolver * restrict self)
{
        if (self) {
                FreeAndClear(self->hessian);
                FreeAndClear(self->gradient);
                FreeAndClear(self);
        }
}

/*
 * Returns the full expected Hessian, with indices corresponding to the
 * full gradient. The matrix is unpacked.
 */
static double * ExpectedHessian(WeightSolver * restrict self)
{
        if (self->hessian) return self->hessian;
        const ModelSpace * restrict space = SolutionModelSpace(self->solution);
        const Model * restrict model = ModelSpaceModel(space);
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return NULL;
        const size_t weightsSize = WeightsSize(model);
        if (!weightsSize) return NULL;
        double * restrict coordinateHessian;
        coordinateHessian = NewExpectedCoordinateHessian(self);
        if (!coordinateHessian) return NULL;
        double * restrict jointHessian;
        jointHessian = NewExpectedCoordinateSpecificityHessian(self);
        if (!jointHessian) return coordinateHessian;
        double * restrict specificityHessian;
        specificityHessian = NewExpectedSpecificityHessian(self);
        if (!specificityHessian) {
                FreeAndClear(jointHessian);
                return coordinateHessian;
        }
        self->hessian = SafeMalloc(self->hessianSize, sizeof(double));
        for (size_t x = 0; x < weightsSize; x++) {
                cblas_dcopy((int)weightsSize, 
                            coordinateHessian + weightsSize * x, 
                            1, 
                            self->hessian + self->gradientSize * x, 
                            1);
                cblas_dcopy((int)classCount,
                            jointHessian + classCount * x,
                            1,
                            (self->hessian 
                             + self->gradientSize * x 
                             + weightsSize),
                            1);
        }
        for (size_t s = 0; s < classCount; s++) {
                cblas_dcopy((int)weightsSize,
                            jointHessian + s,
                            (int)classCount,
                            (self->hessian 
                             + self->gradientSize * (weightsSize + s)),
                            1);
                cblas_dcopy((int)classCount,
                            specificityHessian + classCount * s,
                            1,
                            (self->hessian 
                             + self->gradientSize * (weightsSize + s)
                             + weightsSize),
                            1);
        }
        FreeAndClear(specificityHessian);
        FreeAndClear(jointHessian);
        FreeAndClear(coordinateHessian);
        return self->hessian;
}

/**
 * Returns true if the gradient is negative for a frozen weight.
 */
static bool FrozenGradientIsNegative(const WeightSolver * restrict self,
                                     Solution * restrict solution0)
{
        if (!self) return false;
        const Parameters * restrict parameters;
        parameters = SolutionParameters(self->solution);
        const ModelSpace * restrict space = SolutionModelSpace(solution0);
        const double * restrict weights = Weights(space);
        for (size_t w = 0; w < self->gradientSize; w++)
                if (isless(weights[w], 
                           parameters->minWeightValue + parameters->margin)
#ifdef WINSBERG_LEGACY
                    && isless(self->gradient[w], -1e-3))
#else
                    && isless(self->gradient[w], 0.0))
#endif
                        return true;
        return false;
}

/**
 * Returns the search direction for a line search for coordinates: in this case,
 * the (Moore-Penrose pseudo)-inverse Hessian times the gradient.
 */
static double * NewSearchDirection(WeightSolver * restrict self,
                                   Solution * restrict solution0)
{
        if (!self || !self->gradient) return NULL;
        double * hessian = ExpectedHessian(self);
        if (!hessian) return NULL;
        const Parameters * restrict params = SolutionParameters(solution0);
        const ModelSpace * restrict space = SolutionModelSpace(solution0);
        const double * restrict weights = Weights(space);
        size_t reducedGradientSize = self->gradientSize;
        // The following allocation is possibly too large, but it allows the for
        // loop following to fill the permutation at the same time as it 
        // determines the reduced size.
        size_t * restrict permutation = SafeCalloc(self->gradientSize,
                                                   sizeof(size_t));
        size_t nextIndex = 0;
        double minGradValue = 0.0;
        size_t minGradIndex = self->gradientSize;
        for (size_t w = 0; w < self->gradientSize; w++) {
                if (isless(weights[w], 
                           params->minWeightValue + params->margin)) {
                        reducedGradientSize--;
                        if (isless(self->gradient[w], minGradValue)) {
                                minGradIndex = w;
                                minGradValue = self->gradient[minGradIndex];
                        }
                }
                else permutation[nextIndex++] = w;                
        }
        if (isless(minGradValue, 0.0) 
            // Check the gradient size just in case something went wrong.
            && reducedGradientSize < self->gradientSize) {
                permutation[nextIndex] = minGradIndex;
                reducedGradientSize++;
        }        
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
        FreeAndClear(work);
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
        FreeAndClear(work); 
        FreeAndClear(s);        
        double * restrict searchDirection;
        searchDirection = SafeCalloc(self->gradientSize, sizeof(double));
        for (size_t a = 0; a < reducedGradientSize; a++)
                // LAPACK has overwritten invReducedGradient with the solution.
                searchDirection[permutation[a]] = invReducedGradient[a];
        FreeAndClear(reducedHessian);
        FreeAndClear(invReducedGradient);
        FreeAndClear(permutation);
        return searchDirection;
}

// This routine performs the line search. It follows Winsberg directly.
// N.B.: It may return the same solution as the one passed in.
// WARNING: This is a simple line search that relies on starting with a
//          negative slope. With an insufficient number of iterations, it can
//          return a solution with positive slope.
#ifdef WINSBERG_LEGACY
static Solution * NewLineSearchSolution(WeightSolver * restrict self,
                                        Solution * restrict solution0,
                                        double SSR)
#else
static Solution * NewLineSearchSolution(WeightSolver * restrict self,
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
#ifdef WINSBERG_LEGACY
        // Winsberg curiously only normalises the direction for coordinates
        // but still uses the normalised slope for computing search directions.
        // The following code is just so that we match.
        cblas_dscal((int)self->gradientSize, norm, searchDirection0, 1);
#else
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
        const double * restrict weights = Weights(space0);
        double maxStepSize = INFINITY;
        for (size_t w = 0; w < self->gradientSize; w++) {
                if (isgreaterequal(searchDirection0[w], 0.0)) continue;
                const double thisMax = ((parameters->minWeightValue 
                                         - weights[w])
                                        / searchDirection0[w]);
                maxStepSize = fmin(maxStepSize, thisMax);
        }
        if (islessequal(maxStepSize, 0.0)) return solution0;
        stepSizes[4] = fmin(stepSizes[4], maxStepSize);
        Solution * restrict solution = solution0;
        WeightSolver * restrict solver = NULL;
        double * restrict nextWeights = SafeMalloc(self->gradientSize,
                                                   sizeof(double));
        for (size_t i = 0; i < parameters->maxSearchIterationCount; i++) {
                cblas_dcopy((int)self->gradientSize, 
                            weights, 
                            1, 
                            nextWeights, 
                            1);
                cblas_daxpy((int)self->gradientSize, 
                            stepSizes[4], 
                            searchDirection0, 
                            1, 
                            nextWeights,
                            1);
                Solution * restrict newSolution;
                newSolution = NewSolutionByUpdatingWeights(solution0,
                                                           nextWeights);
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
                DeleteWeightSolver(solver);
                solver = NewWeightSolver(solution);
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
                if (stepSizes[4] 
                    == maxStepSize * parameters->stepIncreaseFactor)
                        break;
                stepSizes[4] = fmin(stepSizes[4], maxStepSize);
        }
        if (parameters->verbosity >= VERY_VERY_VERY_VERBOSE)
                fprintf(stdout, "\n");
        FreeAndClear(nextWeights);
        DeleteWeightSolver(solver);
        FreeAndClear(searchDirection0);
        return solution;
}

Solution * NewWeightSolution(Solution * restrict solution0)
{
        if (!solution0) return NULL;
        const ModelSpace * restrict space = SolutionModelSpace(solution0);
        const Model * restrict model = ModelSpaceModel(space);
        const size_t classCount = ClassCount(model);
        if (classCount < 2) return solution0;
        const Parameters * restrict params = SolutionParameters(solution0);
        WeightSolver * restrict solver = NewWeightSolver(solution0);
        Solution * restrict soln = solution0;
        double SSR = SumOfSquaredModelError(soln);
        if (params->verbosity >= VERY_VERY_VERBOSE)
                fprintf(stdout,
                        "                \n"
                        "                Optimisation of weights\n"
                        "                \n"
                        "                Squared model error after iteration"
                        " 0: %e\n",
                        SSR);
#ifdef WINSBERG_LEGACY
        WeightSolver * restrict solver0 = solver;
        size_t iccount = 0;
#endif
        for (size_t i = 0; i < params->maxWeightIterationCount; i++) {
#ifdef WINSBERG_LEGACY
        Winsberg:
                if (++iccount >= 5) break;
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
                        if (solver != solver0) DeleteWeightSolver(solver);
#else
                        DeleteWeightSolver(solver);
#endif
                        solver = NewWeightSolver(soln);
                        if (isless(relativeImprovement, 
                                   params->weightImprovementFactor)) {
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
                        // We would have hit break if newSln == solution0
                        DeleteSolutionPreservingClassAssignment(newSolution);
                        break;
                }
#ifndef WINSBERG_LEGACY
                SSR = newSSR;
#endif
        }
        if (params->verbosity >= VERY_VERY_VERBOSE) fprintf(stdout, "\n");
#ifdef WINSBERG_LEGACY
        if (solver0 != solver) DeleteWeightSolver(solver0);
#endif
        DeleteWeightSolver(solver);
        Solution * finalSolution = NewSolutionByNormalisingWeights(soln);
        if (soln != solution0)
                DeleteSolutionPreservingClassAssignment(soln);
        return finalSolution;
}
