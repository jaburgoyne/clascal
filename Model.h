// CLASCAL (Model.h)
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
 * @file Model.h
 *
 * Model type and related definitions.
 */

struct _Model;
/**
 * A latent-class model recording the number latent classes, the number of 
 * dimensions, and the type of stimulus specificities.
 */
typedef struct _Model Model;

/**
 * The type of specificities encoded in the model.
 */
enum _SpecificityType {
        /** No specificities are included. */
        NoSpecificities = 0,
        /** Each stimulus has a particular specificity. */ 
        GlobalSpecificities = 1,
        /** Each class has its own specificity per stimulus. */
        ClassSpecificities = 2
};
/**
 * The type of specificities encoded in the model.
 */
typedef enum _SpecificityType SpecificityType;

/**
 * Allocates, initialises, and returns a new Model.
 *
 * @param dimensionCount  the number of Euclidean dimensions (greater than zero)
 * @param classCount  the number of latent classes (greater than zero)
 * @param specificityType  the type of specificties to use
 */
Model * NewModel(size_t, size_t, SpecificityType);

/**
 * Deallocates a Model.
 */
void DeleteModel(Model *);

/** 
 * Returns the number of latent classes in the model. In Winsberg and De Soete
 * 1993, this value is represented by \f$T\f$.
 */
size_t ClassCount(const Model *);

/**
 * Returns the number of dimensions in the model. In Winsberg and De Soete 1993,
 * this value is represented by \f$R\f$.
 */
size_t DimensionCount(const Model *);

/**
 * Returns the type of stimulus specificities in the model.
 */
SpecificityType ModelSpecificityType(const Model *);

/**
 * Returns true if the model has a spatial component (either dimensions
 * or specificities). If false, the model represents the 'null' model for the
 * specified number of classes, i.e., a model with 
 */
bool IsSpatial(const Model *);
