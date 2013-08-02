// CLASCAL (_Solution.h)
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
 * Allocates, initialises, and returns a new Solution. N.B.: The Solution
 * seizes ownership of its ModelSpace and ClassAssignment: they will be 
 * deleted when the Solution is deleted.
 *
 * @param experiment  the experiment used for this solution
 * @param assignment  the class assignment in this solution
 * @param space  the spatial model in this solution
 * @param prior  the a priori distribution over classes
 * @param sigmaSquared  the estimated variance of predictions from this solution
 * @param parameters  the parameters used to generate this solution
 */
Solution * NewSolution(const Experiment *,
                       ClassAssignment *,
                       ModelSpace *,
                       const double *,
                       double,
                       const Parameters *);

/**
 * Returns the estimated number of subjects in each class who have rated each
 * stimulus pair, given the current class assignment and experimental data.
 * Rows are subjects and columns are stimulus pairs. Necessary as an alternative
 * to overall class sizes in case of missing ratings in the data.
 */
const double * PairwiseClassSizes(const Solution *);

/**
 * Returns the weighted mean dissimilarity ratings for each class given 
 * this solution. In Winsberg and De Soete 1993, this matrix is represented by
 * \f$\mathbf{\bar{Y}}\f$. Row indices correspond to latent classes and column
 * indices correspond to stimulus pairs.
 */
const double * ClassDissimilarities(const Solution *);

/**
 * Returns an array of all class residuals for this solution. In Winsberg and 
 * De Soete 1993, this is represented by \f$\bar{y}_{tjk} - \delta_{tjk}\f$. Row 
 * indices are classes and column indices are stimulus pairs.
 */
const double * ClassResiduals(const Solution *);

/**
 * Returns the relative errors: residuals over modelled distances for each 
 * class. Row indices are classes and column indices are stimulus pairs.
 */
const double * RelativeErrors(Solution *);

/**
 * Returns the relative errors weighted by class size.
 */
const double * WeightedRelativeErrors(Solution *);

/**
 * Returns a matrix of values used often in computing the expected Hessians:
 * the class sizes divided by the square of the estimated distances for each
 * class. Row indices correspond to classes and column indices to stimulus
 * pairs.
 */
const double * HessianFactors(Solution *);

/**
 * Allocates, initialises, and returns a new Solution based on a given
 * Solution but with new weights.
 *
 * @param solution  the Solution upon which the new Solution will be based
 * @param weights  the weights used for each class in the Solution, which may be
 *                 NULL in the case of a single class and otherwise has row
 *                 indices corresponding to dimensions and column indices to
 *                 classes (the transpose of the published version)
 */
Solution * NewSolutionByUpdatingWeights(const Solution *, const double *);

/**
 * Allocates, initialises, and returns a new Solution based on a given
 * Solution but with normalised weights.
 *
 * @param solution  the Solution upon which the new Solution will be based
 */
Solution * NewSolutionByNormalisingWeights(const Solution *);

/**
 * Allocates, initialises, and returns a new Solution based on a given
 * Solution but with new coordinates and specificities.
 *
 * @param solution  the Solution upon which the new Solution will be based
 * @param coordinates  the coordinates used in the Solution for each stimulus,
 *                     row indices corresponding to stimuli and column indices
 *                     to dimensions
 * @param specificities  the specificity values (or NULL for models without
 *                       specificities), row indices corresponding to latent
 *                       classes and column indices to stimuli 
 */
Solution * NewSolutionByUpdatingCoordinates(const Solution *,
                                            const double *,
                                            const double *);

/**
 * Deallocates a Solution without deallocating its class assignment (which is 
 * necessary during the M-step).
 */
void DeleteSolutionPreservingClassAssignment(Solution *);

/**
 * Deallocates a Solution without deallocating its model space.
 */
void DeleteSolutionPreservingModelSpace(Solution *);
