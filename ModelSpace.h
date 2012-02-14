// CLASCAL (ModelSpace.h)
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
 *  @file ModelSpace.h
 *
 * Model space type and related definitions.
 */

struct _ModelSpace;
/**
 * A fit model on a given stimulus set.
 */
typedef struct _ModelSpace ModelSpace;

/**
 * Allocates, initialises, and returns a new ModelSpace.
 *
 * @param stimulusSet  the stimulus set embedded in this space
 * @param model  the model used to represent this space
 * @param weights  the weights used for each class in the space, row indices
 *                 corresponding to dimensions and column indices to classes
 *                 (the transpose of the published version); passing NULL will
                   initialise all weights to 1.0
 * @param coordinates  the coordinates used in the space for each stimulus,
 *                     row indices corresponding to stimuli and column indices
 *                     to dimensions; passing NULL will initialise all 
                       coordinates at random (uniform from -1.0 to 1.0)
 * @param specificities  the specificity values, row indices corresponding to 
 *                       latent classes and column indices to stimuli; passing
 *                       NULL will initialise all specificities to 1.0
 *
 * The weights, coordinates, and specificities are copied. It is safe to delete
 * them after the method call. 
 */
ModelSpace * NewModelSpace(const StimulusSet *,
                           const Model *,
                           const double *,
                           const double *,
                           const double *);

/**
 * Deallocates a ModelSpace.
 */
void DeleteModelSpace(ModelSpace *);

/**
 * Returns the stimulus set embedded in this space.
 */
const StimulusSet * ModelSpaceStimulusSet(const ModelSpace *);

/**
 * Returns the model used to parameterise this space.
 */
const Model * ModelSpaceModel(const ModelSpace *);

/**
 * Returns the weights of each class for the respective dimensions (and
 * specificities, where class-based). In Winsberg and De Soete 1993, this matrix
 * is represented by \f$\mathbf{W}\f$. Row indices correspond to dimensions and 
 * column indices to classes. Note that this arrangement is the transpose of the
 * presentation in the paper so as to keep specificity weights together. 
 * If present, the specificity weights follow the dimensions as an additional 
 * row of the matrix.
 */
const double * Weights(const ModelSpace *);

/**
 * Returns the coordinates for each stimulus. In Winsberg and De Soete 1993,
 * this matrix is represented by \f$\mathbf{X}\f$. Row indices correspond to
 * stimuli and column indices to dimensions.
 */
const double * Coordinates(const ModelSpace *);

/**
 * Returns the specificity values, if any. In the case of class-based
 * specificities, the return value is a matrix the row indices of which
 * correspond to latent classes and the column indices of which to stimuli.
 */
const double * Specificities(const ModelSpace *);

/**
 * Returns the number of active parameters in this model.
 */
size_t ParameterCount(const ModelSpace *);

/**
 * Returns the distance between a pair of stimuli for a given latent class. The
 * returned value includes all relevant weights, i.e., it is the expected
 * dissimilarity rating for a subject in this class given the fit model.
 */
double DistanceForPairAndClass(const ModelSpace *,
                               const StimulusPair *,
                               size_t);
