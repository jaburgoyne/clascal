// CLASCAL (Solution.h)
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

/**
 * @file Solution.h
 *
 * Expectation maximiser and related definitions.
 */

struct _Solution;
/**
 * A potential solution (from an iteration of the EM algorithm). 
 */
typedef struct _Solution Solution;

struct _Parameters {
        /**
         * Whether the coordinates are initialised with classical MDS.
         */
        bool useTorgersonStart;
        /**
         * Whether the subjects are pre-classified with k-means clustering.
         */
        bool useKMeansStart;        
        /**
         * Maximum number of iterations of the EM algorithm overall.
         */
        size_t maxEMIterationCount;
        /**
         * Amount of improvement in log-likelihood necessary for convergence in
         * the EM algorithm overall.
         */
        double EMImprovementAmount;        
        /**
         * Maximum number of iterations of the M-step (i.e., the maximum number
         * of times to update the spatial parameters and weights).
         */
        size_t maxMIterationCount;
        /**
         * Factor of improvement in error necessary for overall convergence in
         * the M-step.
         */
        double MImprovementFactor;
        /** 
         * Maximum number of iterations to use within M-steps when updating
         * spatial coordinates and specificities. 
         */
        size_t maxSpatialIterationCount;
        /** 
         * Minimum acceptable value for specificities. 
         */
        double minSpecificityValue;
        /**
         * Factor of improvement in error necessary for convergence in the M-
         * step for coordinates and specificities.
         */
        double spatialImprovementFactor;
        /** 
         * Maximum number of iterations to use within M-steps when updating
         * spatial coordinates and specificities. 
         */
        size_t maxWeightIterationCount;
        /** 
         * Minimum acceptable value for weights. 
         */
        double minWeightValue;
        /**
         * Factor of improvement in error necessary for convergence in the M-
         * step for weights.
         */
        double weightImprovementFactor;
        /** 
         * Maximum number of iterations to use for line search.
         */
        size_t maxSearchIterationCount;
        /** 
         * Step size of the first step during line search. 
         */
        double firstStepSize;
        /**
         * The distance from the boundary at which constrained parameters
         * become frozen.
         */
        double margin;
        /** 
         * Reciprocal of the worst possible condition number allowable for
         * solving the linear-least-squares problem during the line search.
         */
        double invConditionNumber;
        /** 
         * Factor of improvement in slope necessary for convergence in
         * line searches.
         */
        double slopeImprovementFactor;
        /** Factor of increase in slope between line-search iterations. */
        double stepIncreaseFactor;
        /** Factor of decrease in slope between line-search iterations. */
        double stepDecreaseFactor;
};
/**
 * Parameters used for generating a solution.
 */
typedef struct _Parameters Parameters;

#ifdef WINSBERG_LEGACY
static const Parameters DEFAULT_PARAMETERS = {
        false,
        false,
        150,
        1e-3,
        20,
        1e-2,
        5, // In principle, this could be as high as 7, but most tests used 5.
        0.0,
        5e-3, 
        4, // The tests use 5, but SW's code prevents more than 4 from running
        0.0,
        5e-3,
        5, 
        1.0,
        5e-3,
        1e-9, // -log of this should be equal to desired digits of precision
        0.1,
        2.0,
        0.5
};
#else
static const Parameters DEFAULT_PARAMETERS = {
        true,
        true,
        150,
        1e-3,
        20,
        1e-2,
        5, 
        0.0,
        5e-3,
        5,
        0.0,
        5e-3,
        5, 
        1.0,
        5e-3,
        1e-9, // -log of this should be equal to desired digits of precision
        0.1,
        2.0,
        0.5
};
#endif

/**
 * Allocates, initialises, and returns the CLASCAL solution for a given 
 * experiment, model, and parameters.
 */
Solution * NewSolutionForExperimentAndModel(const Experiment *,
                                            const Model *,
                                            const Parameters *);

/**
 * Allocates, initialises, and returns a saved Solution from a file. N.B. The 
 * current implementation leaks memory when loading the parameters.
 *
 * @param file  the file where the Solution was saved
 * @param experiment  the experiment used to compute the solution
 */
Solution * NewSolutionFromFile(FILE *, const Experiment *);

/**
 * Allocates, initialises, and returns a saved Solution from a filename. N.B.
 * The current implementation leaks memory when loading the parameters.
 *
 * @param file  the path of the file where the Solution was saved
 * @param experiment  the experiment used to compute the solution
 */
Solution * NewSolutionFromFilename(const char *, const Experiment *);

/**
 * Deallocates a Solution.
 */
void DeleteSolution(Solution *);

/**
 * Saves a solution to a file.
 */
void SaveSolutionToFile(const Solution *, FILE *);

/**
 * Saves a solution with a given filename.
 */
void SaveSolutionToFilename(const Solution *, const char *);

/**
 * Returns the experiment used to generate a solution.
 */
const Experiment * SolutionExperiment(const Solution *);

/**
 * Returns the class assignment suggested by a solution.
 */
const ClassAssignment * SolutionClassAssignment(const Solution *);

/**
 * Returns the model space suggested by a solution.
 */
const ModelSpace * SolutionModelSpace(const Solution *);

/**
 * Returns the prior probabilities of membership in each class. This is the
 * vector \f$\boldsymbol{\lambda}\f$ in Winsberg and De Soete 1993.
 */
const double * PriorDistribution(const Solution *);

/**
 * Returns the parameters used for computing a solution.
 */
const Parameters * SolutionParameters(const Solution *);

/**
 * Returns the degrees of freedom (data size minus number of parameters) of this 
 * solution.
 */
size_t DegreesOfFreedom(const Solution *);

/**
 * Returns the sum of squared errors of empirical class means versus the
 * predicted model. This is named \f$q^*_1\f$ in Winsberg and De Soete 1993. 
 */
double SumOfSquaredModelError(const Solution *);

/**
 * Returns the estimated variance for a solution.
 */
double EstimatedVariance(const Solution *);

/**
 * Returns the log likelihood of a solution.
 */
double LogLikelihood(const Solution *);

/**
 * Returns Akaike's information criterion for a solution.
 */
double AkaikeCriterion(const Solution *);

/**
 * Returns the Bayesian information criterion for a solution.
 */
double BayesianCriterion(const Solution *);

/**
 * Returns a new experiment with random data generated according to the 
 * parameters in this model.
 */
Experiment * NewMonteCarloExperiment(const Solution *);
