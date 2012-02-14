// CLASCAL (_ModelSpace.h)
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
 * Allocates, initialises, and returns a dummy ModelSpace with artificial
 * distances.
 *
 * @param stimulusSet  the stimulus set embedded in this space
 * @param model  the model used to represent this space
 * @param distances  the artificial distances
 *
 * N.B. The new ModelSpace takes over memory management of the distances array.
 */
ModelSpace * NewDummySpace(const StimulusSet *, const Model *, double *);

/**
 * Returns an array of all distances between stimulus pairs, weighted for each
 * class. Row indices correspond to classes and column indices correspond to
 * stimulus pairs.
 */
const double * ClassDistances(const ModelSpace *);

/**
 * Returns an array of raw difference vectors between the stimulus pairs in the
 * space. Row indices correspond to stimulus pairs and column indices correspond
 * to dimensions.
 */
const double * RawDifferences(const ModelSpace *);

/**
 * Returns the raw differences squared.
 */
const double * SquaredDifferences(const ModelSpace *);

/**
 * Returns the pairwise specificity values (i.e., the sum for each stim. pair).
 */
const double * PairwiseSpecificities(const ModelSpace *);

/**
 * Returns the total number of coordinates (stimulus count times dimension 
 * count).
 */
size_t CoordinatesSize(const ModelSpace *);

/**
 * Returns the total number of specificity values (depending on the specificity
 * type, 0, the stimulus count, or the stimulus count times the class count).
 */
size_t SpecificitiesSize(const ModelSpace *);

/**
 * Returns the total number of pairwise differences for each coordinate: pair 
 * count times dimension count.
 */
size_t DifferencesSize(const ModelSpace *);

/**
 * Returns the total number of modelled dissimilarities: pair count times class 
 * count.
 */
size_t DistancesSize(const ModelSpace *);

/**
 * Allocates, initialises, and returns a new ModelSpace based on a given
 * ModelSpace but with new weights.
 *
 * @param modelSpace  the model space upon which the new space will be based
 * @param weights  the weights used for each class in the space, which may be
 *                 NULL in the case of a single class and otherwise has row
 *                 indices corresponding to dimensions and column indices to
 *                 classes (the transpose of the published version)
 */
ModelSpace * NewModelSpaceByUpdatingWeights(const ModelSpace *,
                                            const double *);

/**
 * Allocates, initialises, and returns a new ModelSpace based on a given
 * ModelSpace but with normalised weights.
 *
 * @param modelSpace  the model space upon which the new space will be based
 */
ModelSpace * NewModelSpaceByNormalisingWeights(const ModelSpace *);

/**
 * Allocates, initialises, and returns a new ModelSpace based on a given
 * ModelSpace but with new coordinates and specificities.
 *
 * @param modelSpace  the model space upon which the new space will be based
 * @param coordinates  the coordinates used in the space for each stimulus,
 *                     row indices corresponding to stimuli and column indices
 *                     to dimensions
 * @param specificities  the specificity values (or NULL for models without
 *                       specificities), row indices corresponding to latent
 *                       classes and column indices to stimuli 
 */
ModelSpace * NewModelSpaceByUpdatingCoordinates(const ModelSpace *,
                                                const double *,
                                                const double *);
