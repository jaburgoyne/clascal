// CLASCAL (Solution.c)
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
#include "SubjectSet.h"
#include "_SubjectSet.h"
#include "Experiment.h"
#include "_Experiment.h"
#include "Model.h"
#include "ClassAssignment.h"
#include "_ClassAssignment.h"
#include "ModelSpace.h"
#include "_ModelSpace.h"
#include "Solution.h"
#include "_Solution.h"
#include "_CoordinateSolver.h"
#include "_WeightSolver.h"

#ifdef WINSBERG_LEGACY
#define WINSBERG_LIKELIHOOD
#endif

typedef __CLPK_integer LPInteger;

struct _Solution {
        const Experiment * restrict experiment;
        ClassAssignment * restrict assignment;
        ModelSpace * restrict space;
        size_t degreesOfFreedom;
        double * prior;
        double estimatedVariance;
        const Parameters * restrict parameters;
        double * restrict classDissimilarities;
        double * restrict classResiduals;
        double sumOfSquaredModelError;
        double * restrict squaredPredictionErrors;
        double * restrict unnormalisedMixtures;
        double * restrict unnormalisedClassMixtures;
        double logLikelihood;
        double akaikeCriterion;
        double bayesianCriterion;
        // For lazy initialisations.
        double * restrict relativeErrors;
        double * restrict weightedRelativeErrors;
        double * restrict hessianFactors;
};

/* 
 * Returns random numbers from the standard normal distribution according to the
 * Marsaglia polar method.
 */
static double RandomNormal()
{
        static bool shouldRegenerate = true;
        static double y;
        static double s = 0.0;
        if (shouldRegenerate) {
                double u;
                double v;
                do {
                        u = (-1.0 
                             + (2.0 
                                * (double)rand() 
                                / (double)RAND_MAX));
                        v = (-1.0 
                             + (2.0 
                                * (double)rand() 
                                / (double)RAND_MAX));
                        s = u * u + v * v;
                } while (s >= 1.0);
                const double multiplier = sqrt(-2.0 * log(s) / s);
                y = multiplier * v;
                shouldRegenerate = false;
                return multiplier * u;
        } else {
                shouldRegenerate = true;
                return y;
        }
}

static double * NewClassDissimilarities(const Solution * restrict self)
{
        const double * restrict classSizes = ClassSizes(self->assignment);
        if (!classSizes) return NULL;
        const double * restrict conditionalDistributions;
        conditionalDistributions = SubjectClassDistributions(self->assignment);
        if (!conditionalDistributions) return NULL;
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        if (!subjectCount) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const Model * restrict model;
        model = ModelSpaceModel(self->space);
        const size_t classCount = ClassCount(model);
        if (!classCount) return NULL;
        const double * restrict subjectDissimilarities;
        subjectDissimilarities = SubjectDissimilarities(self->experiment);
        if (!subjectDissimilarities) return NULL;
        
        const size_t dissimilaritiesSize = DistancesSize(self->space);
        double * restrict classDissimilarities = SafeCalloc(dissimilaritiesSize,
                                                            sizeof(double));
        cblas_dgemm(CblasRowMajor, 
                    CblasTrans, 
                    CblasNoTrans, 
                    (int)classCount, 
                    (int)pairCount, 
                    (int)subjectCount, 
                    1.0, 
                    conditionalDistributions, 
                    (int)classCount, 
                    subjectDissimilarities, 
                    (int)pairCount, 
                    0.0, 
                    classDissimilarities, 
                    (int)pairCount);
        for (size_t i = 0; i < classCount; i++)
                cblas_dscal((int)pairCount, 
                            1.0 / classSizes[i], 
                            classDissimilarities + pairCount * i, 
                            1);
        return classDissimilarities;
}

static double * NewClassResiduals(const Solution * restrict self)
{
        if (!self->classDissimilarities) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!(pairCount)) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!(classCount)) return NULL;
        const double * restrict distances = ClassDistances(self->space);
        if (!distances) return NULL;
        const size_t residualsSize = DistancesSize(self->space);
        double * restrict residuals = SafeMalloc(residualsSize, sizeof(double));
        cblas_dcopy((int)residualsSize, 
                    self->classDissimilarities,
                    1, 
                    residuals,
                    1);
        cblas_daxpy((int)residualsSize, 
                    -1.0, 
                    distances, 
                    1,
                    residuals, 
                    1);
        return residuals;
}

/* N.B.: Class residuals must be initialised before calling this function. */
static double TotalModelError(const Solution * restrict self)
{
        if (!self->classResiduals) return NAN;
        const double * classSizes = ClassSizes(self->assignment);
        if (!classSizes) return NAN;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NAN;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NAN;       
        const size_t residualsSize = DistancesSize(self->space);
        double * restrict squaredResiduals = SafeMalloc(residualsSize, 
                                                        sizeof(double));
        cblas_dcopy((int)residualsSize,
                    self->classResiduals,
                    1,
                    squaredResiduals,
                    1);
        cblas_dtbmv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)residualsSize, 
                    0,
                    self->classResiduals, 
                    1, 
                    squaredResiduals,
                    1); // trick for elementwise multiplication         
        double * restrict accumulator = SafeCalloc(pairCount, sizeof(double));
        cblas_dgemv(CblasRowMajor, 
                    CblasTrans, 
                    (int)classCount, 
                    (int)pairCount, 
                    1.0, 
                    squaredResiduals,
                    (int)pairCount, 
                    classSizes, 
                    1, 
                    0.0, 
                    accumulator, 
                    1);
        free(squaredResiduals);
        // Absolute value should be safe because all values should be positive.
        double sum = cblas_dasum((int)pairCount,
                                 accumulator, 
                                 1);
        free(accumulator);
        return sum;
}

/**
 * Returns the squared errors between subject data and predicted class means.
 * Rows contain subjects and columns classes.
 */
static double * NewSquaredPredictionErrors(const Solution * restrict self)
{
        const double * restrict subjectDissimilarities;
        subjectDissimilarities = SubjectDissimilarities(self->experiment);
        if (!subjectDissimilarities) return NULL;
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        if (!subjectCount) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NULL;
        const double * classDistances = ClassDistances(self->space);
        if (!classDistances) return NULL;
        const size_t distributionsSize = DistributionsSize(self->assignment);
        if (!distributionsSize) return NULL;
        double * restrict residual = SafeMalloc(pairCount, sizeof(double));
        double * restrict squaredPredictionErrors;
        squaredPredictionErrors = SafeMalloc(distributionsSize, sizeof(double));
        for (size_t i = 0; i < subjectCount; i++) {
                for (size_t t = 0; t < classCount; t++) {
                        cblas_dcopy((int)pairCount, 
                                    subjectDissimilarities + pairCount * i, 
                                    1, 
                                    residual, 
                                    1); 
                        cblas_daxpy((int)pairCount, 
                                    -1.0, 
                                    classDistances + pairCount * t, 
                                    1, 
                                    residual, 
                                    1);
                        double dot = cblas_ddot((int)pairCount, 
                                                residual, 
                                                1, 
                                                residual,
                                                1);
                        squaredPredictionErrors[classCount * i + t] = dot;
                }
        }
        free(residual);
        return squaredPredictionErrors;;
}

/** 
 * Returns the probability densities labelled as \f$h\f$ in Winsberg and De 
 * Soete 1993 without the normalisation term. N.B.: The prediction errors must 
 * be initialised before calling as well as the estimated variance. Rows contain
 * subjects and columns contain classes.
 */
static double * NewUnnormalisedMixtures(const Solution * restrict self)
{
        if (!self->squaredPredictionErrors || !self->prior) return NULL;
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        if (!subjectCount) return NULL;        
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NULL;
        const size_t distributionsSize = DistributionsSize(self->assignment);
        if (!distributionsSize) return NULL;
        double * restrict unnormalisedMixtures = SafeMalloc(distributionsSize, 
                                                            sizeof(double));
        cblas_dcopy((int)distributionsSize, 
                    self->squaredPredictionErrors, 
                    1, 
                    unnormalisedMixtures, 
                    1);
        cblas_dscal((int)distributionsSize, 
                    -0.5 / self->estimatedVariance,
                    unnormalisedMixtures, 
                    1);
        for (size_t i = 0; i < distributionsSize; i++)
#ifdef WINSBERG_LIKELIHOOD
                unnormalisedMixtures[i] = (unnormalisedMixtures[i] > 1500
                                           ? 0.0
                                           : exp(unnormalisedMixtures[i]));
#else
                unnormalisedMixtures[i] = exp(unnormalisedMixtures[i]);
#endif
        for (size_t t = 0; t < classCount; t++)                
                cblas_dscal((int)subjectCount,
#ifdef WINSBERG_LIKELIHOOD
                            (self->prior[t] > 10e-30) ? self->prior[t] : 0.0,
#else
                            self->prior[t],
#endif
                            unnormalisedMixtures + t,
                            (int)classCount);
        return unnormalisedMixtures;
}

/**
 * Returns the unnormalised mixtures summed over classes. N.B.: The plain
 * unnormalised mixtures must obviously be initialised first.
 */
static double * NewClassMixtures(const Solution * restrict self)
{
        if (!self->unnormalisedMixtures) return NULL;
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        if (!subjectCount) return NULL;        
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NULL;
        double * restrict unnormalisedClassMixtures;
        unnormalisedClassMixtures = SafeMalloc(subjectCount, sizeof(double));
        for (size_t i = 0; i < subjectCount; i++) {
                // The sum of absolute values should be acceptable because all
                // values should be non-negative.
                const double sum = cblas_dasum((int)classCount, 
                                               (self->unnormalisedMixtures
                                                + classCount * i), 
                                               1);
                unnormalisedClassMixtures[i] = sum;
        }
        return unnormalisedClassMixtures;
}

/* N.B.: The unnormalised class mixtures must be initialised before calling. */
static double LogLikelihoodValue(const Solution * restrict self)
{
        if (!self->unnormalisedClassMixtures) return NAN;
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        if (!subjectCount) return NAN;        
        double logLikelihood = (-0.5 
                                * (double)DataSize(self->experiment)
                                * log(2.0 * M_PI * self->estimatedVariance));
        for (size_t i = 0; i < subjectCount; i++) 
#ifdef WINSBERG_LIKELIHOOD
                logLikelihood += (self->unnormalisedClassMixtures[i] > 10e-30
                                  ? log(self->unnormalisedClassMixtures[i])
                                  : 0.0);
#else
                logLikelihood += log(self->unnormalisedClassMixtures[i]);
#endif
        return logLikelihood;
}

static double * NewNormalisedPrior(const Solution * restrict self,
                                   const double * restrict prior)
{
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        double * restrict newLambda = SafeMalloc(classCount, sizeof(double));
        // N.B. All values must be non-negative.
        cblas_dcopy((int)classCount, prior, 1, newLambda, 1);
        const double sum = cblas_dasum((int)classCount, newLambda, 1);
        cblas_dscal((int)classCount, 1.0 / sum, newLambda, 1);
        return newLambda;
}

static void InitialiseBasicDerivedValues(Solution * restrict self) {
        if (self) {
                const size_t dataSize = DataSize(self->experiment);
                const size_t parameterCount = ParameterCount(self->space);
                if (dataSize < parameterCount)
                        ExitWithError("Model has more parameters than data");
                self->degreesOfFreedom = dataSize - parameterCount;
                self->classDissimilarities = NewClassDissimilarities(self);
                self->classResiduals = NewClassResiduals(self);
                self->sumOfSquaredModelError = TotalModelError(self);
                (self
                 ->squaredPredictionErrors) = NewSquaredPredictionErrors(self);
                self->unnormalisedMixtures = NewUnnormalisedMixtures(self);
                self->unnormalisedClassMixtures = NewClassMixtures(self);
                self->logLikelihood = LogLikelihoodValue(self);
                self->akaikeCriterion = (-2.0 * self->logLikelihood
                                         + 2.0 * ParameterCount(self->space));
                self->bayesianCriterion = (-2.0 * self->logLikelihood
                                           + (log((double)DataSize(self->
                                                                   experiment))
                                              * ParameterCount(self->space)));
        }
}

Solution * NewSolution(const Experiment * restrict experiment,
                       ClassAssignment * restrict assignment,
                       ModelSpace * restrict space,
                       const double * restrict prior,
                       double sigmaSquared,
                       const Parameters * restrict parameters)
{
        if (!assignment 
            || !space
            || ClassAssignmentModel(assignment) != ModelSpaceModel(space)
            || !prior) 
                return NULL;
        Solution * restrict self;
        if ((self = malloc(sizeof(Solution)))) {
                self->experiment = experiment;
                self->assignment = assignment;
                self->space = space;
                self->prior = NewNormalisedPrior(self, prior);
                self->estimatedVariance = sigmaSquared;
                self->parameters = parameters;
                InitialiseBasicDerivedValues(self);
                self->relativeErrors = NULL;
                self->weightedRelativeErrors = NULL;
                self->hessianFactors = NULL;                        
        }
        return self;
}

Solution * NewSolutionByUpdatingWeights(const Solution * restrict self, 
                                        const double * restrict weights)
{
        return NewSolution(self->experiment, 
                           self->assignment, 
                           NewModelSpaceByUpdatingWeights(self->space, weights),
                           self->prior,
                           self->estimatedVariance,
                           self->parameters);
}

Solution * NewSolutionByNormalisingWeights(const Solution * restrict self)
{
        return NewSolution(self->experiment, 
                           self->assignment, 
                           NewModelSpaceByNormalisingWeights(self->space),
                           self->prior,
                           self->estimatedVariance,
                           self->parameters);
}

Solution * 
NewSolutionByUpdatingCoordinates(const Solution * restrict self,
                                 const double * restrict coordinates,
                                 const double * restrict specificities)
{
        return NewSolution(self->experiment, 
                           self->assignment, 
                           NewModelSpaceByUpdatingCoordinates(self->space, 
                                                              coordinates, 
                                                              specificities),
                           self->prior,
                           self->estimatedVariance,
                           self->parameters);
}

static void DeleteSolutionPreservingSpaceAndAssignment(Solution * restrict self)
{
        if (self) {
                FreeAndClear(self->hessianFactors);
                FreeAndClear(self->weightedRelativeErrors);
                FreeAndClear(self->relativeErrors);
                FreeAndClear(self->unnormalisedClassMixtures);
                FreeAndClear(self->unnormalisedMixtures);
                FreeAndClear(self->squaredPredictionErrors);
                FreeAndClear(self->classResiduals);
                FreeAndClear(self->classDissimilarities);
                FreeAndClear(self->prior);
                free(self);
        }
}


void DeleteSolutionPreservingModelSpace(Solution * restrict self)
{
        if (self) {
                DeleteClassAssignment(self->assignment);
                self->assignment = NULL;
                DeleteSolutionPreservingSpaceAndAssignment(self);
        }
}

void DeleteSolutionPreservingClassAssignment(Solution * restrict self)
{
        if (self) {
                DeleteModelSpace(self->space);
                self->space = NULL;
                DeleteSolutionPreservingSpaceAndAssignment(self);
        }
}

void DeleteSolution(Solution * restrict self)
{
        if (self) {
                DeleteClassAssignment(self->assignment);
                self->assignment = NULL;
                DeleteModelSpace(self->space);
                self->space = NULL;
                DeleteSolutionPreservingSpaceAndAssignment(self);
        }
}

const Experiment * SolutionExperiment(const Solution * restrict self)
{
        return self ? self->experiment : NULL;
}

const ClassAssignment * SolutionClassAssignment(const Solution * restrict self)
{
        return self ? self->assignment : NULL;
}

const ModelSpace * SolutionModelSpace(const Solution * restrict self)
{
        return self ? self->space : NULL;
}

size_t DegreesOfFreedom(const Solution * restrict self)
{
        return self ? self->degreesOfFreedom : 0;
}

const double * PriorDistribution(const Solution * restrict self)
{
        return self ? self->prior : NULL;
}

const Parameters * SolutionParameters(const Solution * restrict self)
{
        return self ? self->parameters : &DEFAULT_PARAMETERS;
}

const double * ClassDissimilarities(const Solution * restrict self)
{        
        return self ? self->classDissimilarities : NULL;
}

const double * ClassResiduals(const Solution * restrict self)
{
        return self ? self->classResiduals : NULL;
}

double SumOfSquaredModelError(const Solution * restrict self)
{
        return self ? self->sumOfSquaredModelError : NAN;
}

double EstimatedVariance(const Solution * restrict self)
{
        return self ? self->estimatedVariance : NAN;
}

double LogLikelihood(const Solution * restrict self)
{
        return self ? self->logLikelihood : NAN;
}

double AkaikeCriterion(const Solution * restrict self)
{
        return self ? self->akaikeCriterion : NAN;
}

double BayesianCriterion(const Solution * restrict self)
{
        return self ? self->bayesianCriterion : NAN;
}

/* N.B.: For protected calls only due to lazy initialisation. */
const double * RelativeErrors(Solution * restrict self)
{
        if (!self) return NULL;
        if (self->relativeErrors) return self->relativeErrors;
        if (!self->classResiduals) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        const double * restrict distances = ClassDistances(self->space);
        if (!distances) return NULL;
        if (!classCount) return NULL;               
        const size_t relativeErrorsSize = DistancesSize(self->space);
        self->relativeErrors = SafeMalloc(relativeErrorsSize, sizeof(double));
        cblas_dcopy((int)relativeErrorsSize,
                    self->classResiduals,
                    1,
                    self->relativeErrors,
                    1);
        cblas_dtbsv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)relativeErrorsSize,
                    0,
                    distances, 
                    1, 
                    self->relativeErrors, 
                    1);
        return self->relativeErrors;
}

/* N.B.: For protected calls only due to lazy initialisation. */
const double * WeightedRelativeErrors(Solution * restrict self)
{
        if (!self) return NULL;
        if (self->weightedRelativeErrors) return self->weightedRelativeErrors;
        const double * restrict relativeErrors;
        relativeErrors = RelativeErrors(self);
        if (!relativeErrors) return NULL;
        const double * classSizes = ClassSizes(self->assignment);
        if (!classSizes) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NULL;               
        const size_t relativeErrorsSize = DistancesSize(self->space);
        self->weightedRelativeErrors = SafeMalloc(relativeErrorsSize, 
                                                  sizeof(double));
        cblas_dcopy((int)relativeErrorsSize,
                    relativeErrors,
                    1,
                    self->weightedRelativeErrors,
                    1);
        for (size_t t = 0; t < classCount; t++)
                cblas_dscal((int)pairCount, 
                            classSizes[t], 
                            self->weightedRelativeErrors + pairCount * t, 
                            1);
        return self->weightedRelativeErrors;
}

/* N.B.: For protected calls only due to lazy initialisation. */
const double * HessianFactors(Solution * restrict self)
{
        if (!self) return NULL;
        if (self->hessianFactors) return self->hessianFactors;
        const double * classSizes = ClassSizes(self->assignment);
        if (!classSizes) return NULL;
        const double * restrict distances = ClassDistances(self->space);
        if (!distances) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ModelSpaceStimulusSet(self->space);
        const size_t pairCount = StimulusPairCount(stimulusSet);
        if (!pairCount) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        if (!classCount) return NULL;               
        const size_t distancesSize = DistancesSize(self->space); 
        self->hessianFactors = SafeMalloc(distancesSize, sizeof(double));
        for (size_t m = 0; m < pairCount; m++)
                cblas_dcopy((int)classCount, 
                            classSizes, 
                            1, 
                            self->hessianFactors + m, 
                            (int)pairCount);
        cblas_dtbsv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)distancesSize, 
                    0, 
                    distances, 
                    1, 
                    self->hessianFactors, 
                    1);
        cblas_dtbsv(CblasRowMajor, 
                    CblasUpper, 
                    CblasNoTrans, 
                    CblasNonUnit, 
                    (int)distancesSize, 
                    0, 
                    distances, 
                    1, 
                    self->hessianFactors, 
                    1);
        return self->hessianFactors;
}

static double NewVariance(const Solution * restrict self) {
        if (!self || !self->squaredPredictionErrors) return NAN;
        const double * restrict distributions;
        distributions = SubjectClassDistributions(self->assignment);
        if (!distributions) return NAN;
        const size_t distributionsSize = DistributionsSize(self->assignment);
        if (!distributionsSize) return NAN;
        const double totalSquaredError = cblas_ddot((int)distributionsSize, 
                                                    distributions, 
                                                    1, 
                                                    (self
                                                     ->squaredPredictionErrors), 
                                                    1);
        const double sigmaSquared = (totalSquaredError 
                                     / (double)DataSize(self->experiment));
        return sigmaSquared;
}

/**
 * Returns the expected class assignment given an existing solution (i.e.,
 * the expected \f$z_{it}\f$ values).
 */
static ClassAssignment * 
NewExpectedAssignment(const Solution * restrict self)
{
        const SubjectSet * restrict subjectSet;
        subjectSet = ClassAssignmentSubjectSet(self->assignment);
        const Model * restrict model;
        model = ClassAssignmentModel(self->assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        const size_t classCount = ClassCount(model);
        if (subjectCount == classCount) {
                const double * restrict distributions;
                distributions = SubjectClassDistributions(self->assignment);
                return NewClassAssignment(subjectSet, model, distributions);
        }
        if (!self->unnormalisedMixtures) return NULL;
        const size_t distributionsSize = DistributionsSize(self
                                                           ->assignment);
        double * restrict newDistributions = SafeMalloc(distributionsSize,
                                                        sizeof(double));
        cblas_dcopy((int)distributionsSize, 
                    self->unnormalisedMixtures, 
                    1, 
                    newDistributions, 
                    1);
#ifdef WINSBERG_LIKELIHOOD
        for (size_t i = 0; i < subjectCount; i++) {
                for (size_t t = 0; t < classCount; t++) {
                        const size_t it = classCount * i + t;
                        if (islessequal(newDistributions[it], 
                                        10e-60 * self->prior[t]))
                                newDistributions[it] = 0.0;
                }
        }
        double * restrict accumulator;
        accumulator = SafeCalloc(subjectCount, sizeof(double));
        for (size_t t = 0; t < classCount; t++)
                cblas_daxpy((int)subjectCount, 
                            1.0, 
                            newDistributions + t, 
                            (int)classCount, 
                            accumulator,
                            1);
        for (size_t i = 0; i < subjectCount; i++) 
                if (accumulator[i] == 0.0)
                        cblas_dcopy((int)classCount, 
                                    self->prior, 
                                    1, 
                                    newDistributions + classCount * i, 
                                    1);
        free(accumulator);
#endif
        // Normalisation is taken care of within NewClassAssignment.
        ClassAssignment * restrict newAssignment;
        newAssignment = NewClassAssignment(subjectSet, model, newDistributions);
        free(newDistributions);
        return newAssignment;
}

/* 
 * Returns the 'null' (non-spatial) solution based on the starting parameters.
 */
static Solution * NewDummySolution(const Experiment * restrict experiment,
                                   const Model * restrict model,
                                   const Parameters * restrict parameters)
{
        if (!experiment || !model || !parameters) return NULL;
        const size_t classCount = ClassCount(model);
        const SubjectSet * restrict subjectSet;
        subjectSet = ExperimentSubjectSet(experiment);
        const size_t subjectCount = SubjectCount(subjectSet);    
        if (classCount > subjectCount) 
                ExitWithError("Number of classes exceeds number of subjects");
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ExperimentStimulusSet(experiment);
        const size_t stimulusPairCount = StimulusPairCount(stimulusSet);
        const size_t distancesSize = SizeProduct(stimulusPairCount, classCount);
        const size_t distributionsSize = SizeProduct(subjectCount,
                                                     classCount);
        const double * data = SubjectDissimilarities(experiment);
        double * restrict prior = SafeMalloc(classCount, sizeof(double));
        double * restrict initialDistances;
        initialDistances = SafeMalloc(distancesSize, sizeof(double));
        double * restrict initialDistributions;
        if (classCount == subjectCount) {
                // This is INDSCAL, for which obviously the starting points will
                // be the original data points.
                for (size_t t = 0; t < classCount; t++)
                        // Normalisation will happen later.
                        prior[t] = 1.0;
                for (size_t d = 0; d < distancesSize; d++)
                        initialDistances[d] = data[d];
                initialDistributions = SafeCalloc(distributionsSize,
                                                  sizeof(double));
                for (size_t i = 0; i < subjectCount; i++)
                        initialDistributions[classCount * i + i] = 1.0;
        } else if (parameters->useKMeansStart && classCount > 1) {
                const size_t subjectPairCount = SubjectPairCount(subjectSet);
                const SubjectPair * subjectPairs = SubjectPairs(subjectSet);
                const double * squaredDistances;
                squaredDistances = ExperimentSquaredDistances(experiment);
                // Find maximally-dissimilar seeds for each class. We use prior
                // to store the initial seeds to save a malloc call.
                size_t * sizes = SafeMalloc(classCount, sizeof(size_t));
                double maxSquaredDistance = 0.0;
                for (size_t l = 0; l < subjectPairCount; l++) {
                        // N.B. isgreater() will (correctly) ignore NaNs.
                        if (isgreater(squaredDistances[l],
                                      maxSquaredDistance)) {
                                maxSquaredDistance = squaredDistances[l];
                                sizes[0] = subjectPairs[l].i1;
                                sizes[1] = subjectPairs[l].i2;
                        }
                }
                if (maxSquaredDistance == 0.0)
                        ExitWithError("Too many subjects are indistinguisable");
                for (size_t t = 2; t < classCount; t++) {
                        double maxSumOfSquares = 0.0;
                        for (size_t i = 0; i < subjectCount; i++) {
                                bool isSeed = false;
                                for (size_t seed = 0; seed < t; seed++) {
                                        if (i == sizes[seed]) {
                                                isSeed = true;
                                                break;
                                        }                                        
                                }
                                if (isSeed) continue;
                                double sumOfSquares = 0.0;
                                size_t compareCount = 0;
                                for (size_t seed = 0; seed < t; seed++) {                                                
                                        SubjectPair pair = { i, sizes[seed] };
                                        size_t l;
                                        l = PairNumberForSubjectPair(subjectSet, 
                                                                     &pair);
                                        if (!isnan(squaredDistances[l])) {
                                                sumOfSquares
                                                        += squaredDistances[l];
                                                compareCount++;
                                        }
                                }
                                if (compareCount > 0)
                                        sumOfSquares /= (double)compareCount++;
                                else sumOfSquares = NAN;
                                if (isgreater(sumOfSquares, maxSumOfSquares)) {
                                        maxSumOfSquares = sumOfSquares;
                                        sizes[t] = i;
                                }
                        }
                        if (maxSumOfSquares == 0.0)
                                ExitWithError("Too many subjects"
                                              " are indistinguishable");
                }
                for (size_t t = 0; t < classCount; t++)
                        cblas_dcopy((int)stimulusPairCount, 
                                    data + stimulusPairCount * sizes[t], 
                                    1, 
                                    initialDistances + stimulusPairCount * t, 
                                    1);
                // Cluster
                size_t * restrict clustering;
                clustering = SafeCalloc(subjectCount, sizeof(size_t));
                double * restrict accumulator;
                accumulator = SafeMalloc(stimulusPairCount, sizeof(double));
                bool hasChanged = true;
                while (hasChanged) {
                        hasChanged = false;
                        // Assign subjects to clusters and record counts.
                        for (size_t t = 0; t < classCount; t++) sizes[t] = 0;
                        for (size_t i = 0; i < subjectCount; i++) {
                                double minSquaredDistance = INFINITY;
                                size_t cluster = clustering[i];
                                for (size_t t = 0; t < classCount; t++) {
                                        size_t ratingCount = 0;
                                        cblas_dcopy((int)stimulusPairCount,
                                                    (data + 
                                                     stimulusPairCount * i), 
                                                    1, 
                                                    accumulator, 
                                                    1);
                                        cblas_daxpy((int)stimulusPairCount, 
                                                    -1.0, 
                                                    (initialDistances 
                                                     + stimulusPairCount * t), 
                                                    1, 
                                                    accumulator, 
                                                    1);
                                        for (size_t m = 0;
                                             m < stimulusPairCount;
                                             m++) {
                                                if (isnan(accumulator[m]))
                                                        accumulator[m] = 0.0;
                                                else ratingCount++;
                                        }
                                        const double d
                                         = (ratingCount > 0
                                            ? cblas_ddot((int)stimulusPairCount,
                                                         accumulator,
                                                         1,
                                                         accumulator,
                                                         1)
                                            : NAN);
                                        if (isless(d, minSquaredDistance)) {
                                                minSquaredDistance = d;
                                                cluster = t;
                                        }
                                }
                                if (cluster != clustering[i]) {
                                        clustering[i] = cluster;
                                        hasChanged = true;
                                }
                                sizes[clustering[i]]++;
                        }
                        // Recompute means if necessary
                        if (hasChanged) {
                                size_t * restrict ratingCounts;
                                ratingCounts = SafeCalloc(distancesSize,
                                                          sizeof(size_t));
                                for (size_t m = 0; m < distancesSize; m++)
                                        initialDistances[m] = 0.0;
                                for (size_t i = 0; i < subjectCount; i++) {
                                        for (size_t m = 0;
                                             m < stimulusPairCount;
                                             m++) {
                                                const double d
                                                    = data[stimulusPairCount * i
                                                           + m];
                                                const size_t k
                                                    = ((stimulusPairCount
                                                        * clustering[i])
                                                       + m);
                                                if (!isnan(d)) {
                                                    initialDistances[k] += d;
                                                    ratingCounts[k]++;
                                                }
                                        }
                                }
                                for (size_t m = 0; m < distancesSize; m++)
                                        initialDistances[m]
                                                = (ratingCounts[m] > 0
                                                   ? (initialDistances[m]
                                                      /= (double)ratingCounts[m])
                                                   : NAN);
                                free(ratingCounts);
                        }
                }
                free(accumulator);
                initialDistributions = SafeCalloc(distributionsSize,
                                                  sizeof(double));
                for (size_t i = 0; i < subjectCount; i++)
                        initialDistributions[classCount * i 
                                             + clustering[i]] = 1.0;
                free(clustering);
                for (size_t t = 0; t < classCount; t++)
                        // Normalisation will happen later.
                        prior[t] = (double)sizes[t];
                free(sizes);
        } else {
                for (size_t t = 0; t < classCount; t++)
                        // Normalisation will happen later.
                        prior[t] = (double)rand();
                double minD = INFINITY;
                double maxD = -INFINITY;
                for (size_t d = 0; d < DissimilaritiesSize(experiment); d++) {
                        // N.B. fmin() and fmax() handle NaNs properly.
                        minD = fmin(minD, data[d]);
                        maxD = fmax(maxD, data[d]);
                }
                const double range = maxD - minD;
                for (size_t d = 0; d < distancesSize; d++) {
                        initialDistances[d] = (minD 
                                               + (range
                                                  * (double)rand()
                                                  / (double)RAND_MAX));
                }
                initialDistributions = NULL;
        }
        ClassAssignment * restrict dummyAssignment;
        dummyAssignment = NewClassAssignment(subjectSet, 
                                             model, 
                                             initialDistributions);
        free(initialDistributions);
        ModelSpace * restrict dummySpace;
        dummySpace = NewDummySpace(stimulusSet, model, initialDistances);
        // initialDistances is now seized by dummySpace
        Solution * restrict tempSolution;
        tempSolution = NewSolution(experiment, 
                                   dummyAssignment, 
                                   dummySpace, 
                                   prior, 
                                   0.0, // dummy variance to be replace later, 
                                   parameters);
        free(prior);
        const double initialVariance = NewVariance(tempSolution);
        Solution * restrict dummySolution = NULL;
        if ((dummySolution = malloc(sizeof(Solution)))) {
                memcpy(dummySolution, tempSolution, sizeof(Solution));
                dummySolution->prior = SafeMalloc(classCount, sizeof(double));
                cblas_dcopy((int)classCount, 
                            tempSolution->prior, 
                            1, 
                            dummySolution->prior, 
                            1);
                dummySolution->estimatedVariance = initialVariance;
                InitialiseBasicDerivedValues(dummySolution);
        }
        DeleteSolutionPreservingSpaceAndAssignment(tempSolution);
        return dummySolution;
}

/**
 * Returns a new solution with an optimised class assignment (i.e., updates the
 * \f$z_{it}\f$ values). The new model space will be the same as that of the 
 * solution passed as an argument. Be careful to avoid double-freeing these 
 * addresses.
 */
static Solution * NewExpectedSolution(const Solution * restrict self)
{
        return NewSolution(self->experiment, 
                           NewExpectedAssignment(self), 
                           self->space,
                           self->prior,
                           self->estimatedVariance,
                           self->parameters);
}

/* Returns an initial solution based on the starting parameters. */
static Solution * NewInitialSolution(const Experiment * restrict experiment,
                                     const Model * restrict model,
                                     const Parameters * restrict parameters)
{
        if (!experiment || !model || !parameters) return NULL;
        const size_t classCount = ClassCount(model);
        const SubjectSet * restrict subjectSet;
        subjectSet = ExperimentSubjectSet(experiment);
        const size_t subjectCount = SubjectCount(subjectSet);    
        if (classCount > subjectCount) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ExperimentStimulusSet(experiment);
        const size_t stimulusPairCount = StimulusPairCount(stimulusSet);
        const double * data = SubjectDissimilarities(experiment);
        Solution * restrict dummySolution;
        dummySolution = NewDummySolution(experiment, model, parameters);
        // We can circumvent most of the function for non-spatial models.
        if (!IsSpatial(model)) {
                Solution * initialSolution = NewExpectedSolution(dummySolution);
                DeleteSolutionPreservingModelSpace(dummySolution);
                return initialSolution;
        }
        ClassAssignment * restrict initialAssignment;
        initialAssignment = NewExpectedAssignment(dummySolution);
        ModelSpace * restrict initialSpace = NULL;
        if (parameters->useTorgersonStart) {
                const size_t dimensionCount = DimensionCount(model);
                const size_t stimulusCount = StimulusCount(stimulusSet);
                if (dimensionCount > stimulusCount)
                        ExitWithError("Too many dimensions to start using"
                                      " Torgerson's algorithm.");
                // Compute the average squared dissimilarity for each pair.
                double * restrict means;
                means = SafeCalloc(stimulusPairCount, sizeof(double));
#ifdef WINSBERG_LEGACY
                // Winsberg treats zeros specially with more than one class.
                if (subjectCount == 1) {
                        cblas_dcopy((int)stimulusPairCount, 
                                    data, 
                                    1, 
                                    means, 
                                    1);
                        cblas_dtbmv(CblasRowMajor, 
                                    CblasUpper, 
                                    CblasNoTrans, 
                                    CblasNonUnit, 
                                    (int)stimulusPairCount, 
                                    0, 
                                    data, 
                                    1, 
                                    means, 
                                    1); // trick for elementwise multiplication
                } else {
                        double grandMean = 0.0;
                        size_t grandCount = 0;
                        size_t * counts;
                        counts = SafeCalloc(stimulusPairCount, sizeof(size_t));
                        for (size_t m = 0; m < stimulusPairCount; m++) {
                                for (size_t i = 0; i < subjectCount; i++) {
                                        const double d = data[(stimulusPairCount 
                                                               * i)
                                                              + m];
                                        if (isgreater(d, 0.0)) {
                                                means[m] += d * d;
                                                counts[m]++;
                                        }
                                }
                                grandMean += means[m];
                                grandCount += counts[m];
                        }
                        grandMean /= (double)grandCount;
                        for (size_t m = 0; m < stimulusPairCount; m++)
                                means[m] = (counts[m]
                                            ? means[m] / (double)counts[m]
                                            : grandMean);
                        free(counts);
                }
#else
                double * restrict accumulator;
                accumulator = SafeMalloc(stimulusPairCount, sizeof(double));
                for (size_t i = 0; i < subjectCount; i++) {
                        cblas_dcopy((int)stimulusPairCount, 
                                    data + stimulusPairCount * i, 
                                         1, 
                                         accumulator, 
                                         1);
                        cblas_dtbmv(CblasRowMajor, 
                                    CblasUpper, 
                                    CblasNoTrans, 
                                    CblasNonUnit, 
                                    (int)stimulusPairCount, 
                                    0, 
                                    data + stimulusPairCount * i, 
                                    1, 
                                    accumulator, 
                                    1); // trick for elementwise multiplication
                        cblas_daxpy((int)stimulusPairCount, 
                                    1.0 / (double)subjectCount, 
                                    accumulator, 
                                    1, 
                                    means, 
                                    1);
                }
                free(accumulator);
#endif
                // Derive the Gram matrix.
                const size_t gramSize = SizeProduct(stimulusCount,
                                                    stimulusCount);
                double * gram = SafeCalloc(gramSize, sizeof(double));
                const StimulusPair * stimulusPairs = StimulusPairs(stimulusSet);
                for (size_t m = 0; m < stimulusPairCount; m++) {
                        StimulusPair pair = stimulusPairs[m];
                        // We technically only need to store either the upper or
                        // the lower triangle, but averaging rows is easier
                        // when we store the entire thing.
                        gram[stimulusCount * pair.j + pair.k]
                        = gram[stimulusCount * pair.k + pair.j]
                        = means[m];
                }
                for (size_t j = 0; j < stimulusCount; j++)
                        means[j] = cblas_dasum((int)stimulusCount, 
                                               gram + stimulusCount * j, 
                                               1) / (double)stimulusCount;
                for (size_t j = 0; j < stimulusCount; j++)
                        cblas_daxpy((int)stimulusCount, 
                                    -1.0, 
                                    means, 
                                    1, 
                                    gram + stimulusCount * j, 
                                    1);
                for (size_t j = 0; j < stimulusCount; j++)
                        cblas_daxpy((int)stimulusCount, 
                                    -1.0, 
                                    means, 
                                    1, 
                                    gram + j, 
                                    (int)stimulusCount);
                const double grandMean = cblas_dasum((int)stimulusCount, 
                                                     means, 
                                                     1) / (double)stimulusCount;
                for (size_t m = 0; m < gramSize; m++) gram[m] += grandMean;
                cblas_dscal((int)gramSize, -0.5, gram, 1);
                free(means);
                // Run MDS.
                // LAPACK is all call-by-reference.
                const size_t coordinatesSize = SizeProduct(stimulusCount,
                                                           dimensionCount);
                char jobz = 'V';
                char range = 'I';
                char uplo = 'L';
                LPInteger n = (LPInteger)stimulusCount;
                LPInteger lda = (LPInteger)stimulusCount;
                double vl = 0.0;
                double vu = 0.0;
                LPInteger il = (LPInteger)(stimulusCount - dimensionCount + 1);
                LPInteger iu = (LPInteger)stimulusCount;
                char cmach = 'S';
                double abstol = dlamch_(&cmach);
                LPInteger m; // dimensionCount
                double * w = SafeMalloc(stimulusCount, sizeof(double));
                double * z = SafeMalloc(coordinatesSize, sizeof(double));
                LPInteger ldz = (LPInteger)stimulusCount;
                LPInteger * isuppz = SafeMalloc(SizeProduct(2, dimensionCount),
                                                sizeof(LPInteger));
                double * work = SafeMalloc(1, sizeof(double));
                LPInteger lwork = -1;
                LPInteger * iwork = SafeMalloc(1, sizeof(LPInteger));
                LPInteger liwork = -1;
                LPInteger info;
                // This call determines work size only.
                dsyevr_(&jobz, 
                        &range, 
                        &uplo, 
                        &n, 
                        gram, 
                        &lda, 
                        &vl, 
                        &vu, 
                        &il, 
                        &iu, 
                        &abstol, 
                        &m, 
                        w, 
                        z, 
                        &ldz, 
                        isuppz, 
                        work, 
                        &lwork, 
                        iwork, 
                        &liwork, 
                        &info);
                lwork = (LPInteger)lround(work[0]);
                free(work);
                work = SafeMalloc((size_t)lwork, sizeof(double));
                liwork = iwork[0];
                free(iwork);
                iwork = SafeMalloc((size_t)liwork, sizeof(LPInteger));
                // This call is the real thing.
                dsyevr_(&jobz, 
                        &range, 
                        &uplo, 
                        &n, 
                        gram, 
                        &lda, 
                        &vl, 
                        &vu, 
                        &il, 
                        &iu, 
                        &abstol, 
                        &m, 
                        w, 
                        z, 
                        &ldz, 
                        isuppz, 
                        work, 
                        &lwork, 
                        iwork, 
                        &liwork, 
                        &info);
                free(iwork);
                free(work);
                free(isuppz);
                free(gram);
                double * restrict initialCoords;
                initialCoords = SafeMalloc(coordinatesSize, sizeof(double));
                for (size_t r = 0; r < dimensionCount; r++) {
                        cblas_dscal((int)stimulusCount, 
                                    sqrt(w[r]), 
                                    z + stimulusCount * r, 
                                    1);
                        cblas_dcopy((int)stimulusCount,
                                    z + stimulusCount * r,
                                    1,
                                    initialCoords + (dimensionCount - r - 1),
                                    (int)dimensionCount);
                }
                free(z);
                free(w);
                initialSpace = NewModelSpace(stimulusSet, 
                                             model, 
                                             NULL, 
                                             initialCoords, 
                                             NULL);
                free(initialCoords);
        } else {
                initialSpace = NewModelSpace(stimulusSet, 
                                             model, 
                                             NULL, 
                                             NULL,
                                             NULL);                
        }
        Solution * restrict initialSolution;
        initialSolution = NewSolution(experiment, 
                                      initialAssignment, 
                                      initialSpace, 
                                      dummySolution->prior, 
                                      dummySolution->estimatedVariance, 
                                      parameters);
        DeleteSolution(dummySolution);
        return initialSolution;
}

/**
 * Returns a new solution with optimised weights, coordinates, and 
 * specificities. The new class assignment will be the same as the class
 * assignment of the solution passed as an argument. Be careful to avoid 
 * double-freeing these addresses.
 */
static Solution * NewMaximisedSolution(Solution * restrict self)
{
        if (!self) return NULL;
        const Parameters * parameters = SolutionParameters(self);
        const Model * restrict model = ModelSpaceModel(self->space);
        ModelSpace * restrict newSpace;
        if (!IsSpatial(model)) {
                const size_t dissimilaritiesSize = DistancesSize(self->space);
                double * dissimilarities = SafeMalloc(dissimilaritiesSize, 
                                                      sizeof(double));
                cblas_dcopy((int)dissimilaritiesSize, 
                            self->classDissimilarities, 
                            1, 
                            dissimilarities, 
                            1);
                newSpace = NewDummySpace(ModelSpaceStimulusSet(self->space), 
                                         model, 
                                         dissimilarities);
        } else {
                const double * coordinates = Coordinates(self->space);
                const double * specificities = Specificities(self->space);
                newSpace = NewModelSpaceByUpdatingCoordinates(self->space, 
                                                              coordinates, 
                                                              specificities);
        }
        Solution * restrict solution = NewSolution(self->experiment, 
                                                   self->assignment, 
                                                   newSpace,
                                                   self->prior,
                                                   self->estimatedVariance,
                                                   parameters);
        // N.B.: Winsberg uses -1e10 rather than the current error.
        double SSR = SumOfSquaredModelError(solution);
        if (parameters->verbosity >= VERY_VERBOSE)
                fprintf(stdout, 
                        "        \n"
                        "        M-step\n"
                        "        \n"
                        "        Squared model error after iteration 0: %e\n",
                        SSR);
        for (size_t i = 0; i < parameters->maxMIterationCount; i++) {
                // It is inelegant to break here, but it seems worse to 
                // re-organise the whole function with a big if statement.
                if (!IsSpatial(model)) break;
                Solution * restrict spaceSolution;
                spaceSolution = NewSpatialSolution(solution);
                Solution * restrict newSolution;
                newSolution = NewWeightSolution(spaceSolution);
                if (newSolution != spaceSolution && spaceSolution != solution)
                        DeleteSolutionPreservingClassAssignment(spaceSolution);
                double newSSR = SumOfSquaredModelError(newSolution);
                if (parameters->verbosity >= VERY_VERBOSE)
                        fprintf(stdout, 
                                "        Squared model error after iteration"
                                " %zu: %e\n",
                                SizeSum(i, 1),
                                newSSR);
                double relativeImprovement = (SSR - newSSR) / SSR;
                if (newSolution != solution)
                        DeleteSolutionPreservingClassAssignment(solution);
                solution = newSolution;
                SSR = newSSR;
#ifdef WINSBERG_LEGACY
                // Winsberg uses a fake log likelihood at this step that has the
                // wrong sign. The following code replicates the effect.
                if (islessequal(relativeImprovement, 0.0)
                    && isgreater(relativeImprovement, 
                                 parameters->MImprovementFactor))
#else
                if (isgreaterequal(relativeImprovement, 0.0)
                    && isless(relativeImprovement, 
                              parameters->MImprovementFactor))
#endif
                        break; // convergence
        }
        if (parameters->verbosity >= VERY_VERBOSE) fprintf(stdout, "\n");
        const double newVariance = NewVariance(solution);
        const double * restrict newPrior;
        // solution->assignment should be the same as solution0->assignment
        newPrior = UnconditionalClassDistribution(solution->assignment);
        Solution * restrict finalSolution = NewSolution(solution->experiment, 
                                                        solution->assignment, 
                                                        solution->space, 
                                                        newPrior, 
                                                        newVariance, 
                                                        parameters);        
        DeleteSolutionPreservingSpaceAndAssignment(solution);
        return finalSolution;
}

/**
 * Returns the optimal solution starting from the given solution using the EM
 * algorithm. Returns NULL if the algorithm fails to converge.
 */
static Solution * NewOptimalSolution(Solution * restrict self)
{
        if (!self) return NULL;
        const size_t classCount = ClassCount(ModelSpaceModel(self->space));
        const Parameters * parameters = SolutionParameters(self);
        double * restrict logLikelihoods;
        logLikelihoods = SafeMalloc(SizeSum(parameters->maxEMIterationCount, 1),
                                    sizeof(double));
        logLikelihoods[0] = self->logLikelihood;
        if (parameters->verbosity >= VERBOSE)
                fprintf(stdout,
                        "Log likelihood after iteration 0: %e\n",
                        logLikelihoods[0]);
        Solution * restrict solution = self;
        for (size_t i = 1; i <= parameters->maxEMIterationCount; i++) {
                Solution * restrict maxSolution;
                maxSolution = NewMaximisedSolution(solution);
                Solution * restrict newSolution;
                if (classCount > 1) {
                        newSolution = NewExpectedSolution(maxSolution);
                        DeleteSolutionPreservingSpaceAndAssignment(maxSolution);                        
                } else {
                        newSolution = maxSolution;
                }
                if (isnan(newSolution->logLikelihood))
                        ExitWithError("Unstable log likelihood");
                logLikelihoods[i] = newSolution->logLikelihood;
                if (parameters->verbosity >= VERBOSE)
                        fprintf(stdout,
                                "Log likelihood after iteration %zu: %e\n",
                                i,
                                logLikelihoods[i]);
                if (solution != self) {
                        if (classCount > 1) DeleteSolution(solution);
                        else DeleteSolutionPreservingClassAssignment(solution);
                }
                solution = newSolution;
                // Check for convergence (and ensure we reach at least i = 2).
                if ((i > 1
                     && isless(logLikelihoods[i] - logLikelihoods[i-1], 0.0))
                    || (i >= 5 
                        && isless(logLikelihoods[i] - logLikelihoods[i-5],
                                  parameters->EMImprovementAmount)))
                        break; // convergence
        }
        free(logLikelihoods);
        return solution;
}

Solution * 
NewSolutionForExperimentAndModel(const Experiment * restrict experiment,
                                 const Model * restrict model,
                                 const Parameters * restrict parameters)
{
        if (!experiment || !model) return NULL;
        if (ClassCount(model) > SubjectCount(ExperimentSubjectSet(experiment)))
                ExitWithError("Model has more classes than there are subjects");
        const Parameters * guaranteedParameters;
        guaranteedParameters = parameters ? parameters : &DEFAULT_PARAMETERS;
        Solution * restrict initialSolution;
        initialSolution = NewInitialSolution(experiment, 
                                             model, 
                                             guaranteedParameters);
        Solution * restrict solution = NewOptimalSolution(initialSolution);
        if (ClassCount(model) > 1) DeleteSolution(initialSolution);
        else DeleteSolutionPreservingClassAssignment(initialSolution);
        return solution;
}

Experiment * NewMonteCarloExperiment(const Solution * restrict self)
{
        if (!self) return NULL;
        const StimulusSet * restrict stimulusSet;
        stimulusSet = ExperimentStimulusSet(self->experiment);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        char * * restrict blankStimulusNames;
        blankStimulusNames = SafeCalloc(stimulusCount, sizeof(char *));
        StimulusSet * restrict blankStimulusSet;
        blankStimulusSet = NewStimulusSet(stimulusCount, blankStimulusNames);
        const SubjectSet * restrict subjectSet;
        subjectSet = ExperimentSubjectSet(self->experiment);
        const size_t subjectCount = SubjectCount(subjectSet);
        char * * restrict blankSubjectNames;
        blankSubjectNames = SafeCalloc(subjectCount, sizeof(char *));
        SubjectSet * restrict blankSubjectSet;
        blankSubjectSet = NewSubjectSet(subjectCount, blankSubjectNames);
        const size_t stimulusPairCount = StimulusPairCount(stimulusSet);
        size_t classCount = ClassCount(ModelSpaceModel(self->space));
        const double * distances = ClassDistances(self->space);
        double * restrict cumulativePrior;
        cumulativePrior = SafeMalloc(SizeSum(classCount, 1), sizeof(double));
        cumulativePrior[0] = 0.0;
        const double estimatedDeviation = sqrt(self->estimatedVariance);
        for (size_t t = 1; t < classCount; t++)
                cumulativePrior[t] = self->prior[t-1] + cumulativePrior[t-1];
        cumulativePrior[classCount] = 1.0;
        const size_t dissimilaritiesSize = DissimilaritiesSize(self
                                                               ->experiment);
        double * restrict generatedData;
        generatedData = SafeMalloc(dissimilaritiesSize, sizeof(double));
        for (size_t i = 0; i < subjectCount; i++){
                double randomValue = (double)rand() / (double)RAND_MAX;
                size_t class = classCount;
                for (size_t t = 0; t < classCount; t++) {
                        if (randomValue <= cumulativePrior[t+1]) {
                                class = t;
                                break;
                        }
                }
                if (class == classCount)
                        ExitWithError("Unable to generate a random solution");
                for (size_t m = 0; m < stimulusPairCount; m++) {
                        generatedData[stimulusPairCount * i + m]
                        = (distances[stimulusPairCount * class + m]
                           + estimatedDeviation * RandomNormal());
                }
        }        
        free(cumulativePrior);
        return NewExperiment(NULL, 
                             blankStimulusSet, 
                             blankSubjectSet,
                             generatedData);
}
