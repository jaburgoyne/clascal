// CLASCAL (ClassAssignment.h)
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
 * @file ClassAssignment.h
 * 
 * Class assignment type and related definitions.
 */

struct _ClassAssignment;
/**
 * A distribution of subjects from a subject set among classes.
 */
typedef struct _ClassAssignment ClassAssignment;

/**
 * Allocates, initialises, and returns a new ClassAssignment. The distributions
 * will be normalised to sum to one, but all values passed in must be non-
 * negative.
 *
 * @param subjectSet  the subject set
 * @param model  the model
 * @param distributions  the class distributions over classes for each subject,
 *                       row indices corresponding to subjects and column
 *                       indices to classes; passing NULL divides all subjects
 *                       evenly
 */
ClassAssignment * NewClassAssignment(const SubjectSet *,
                                     const Model *,
                                     const double *);

/**
 * Deallocates a ClassAssignment.
 */
void DeleteClassAssignment(ClassAssignment *);

/**
 * Returns the subject set for this class assignment.
 */
const SubjectSet * ClassAssignmentSubjectSet(const ClassAssignment *);

/**
 * Returns the model used for assigning classes.
 */
const Model * ClassAssignmentModel(const ClassAssignment *);

/**
 * Returns the distribution functions over classes for each subject. In Winsberg
 * and De Soete 1993, these values are represented by 
 * \f$\hat{z}_{it} = h_{it}(\mathbf{\hat{X}}, \mathbf{\hat{W}}, \hat{\sigma}^2, 
 * \boldsymbol{\hat{\lambda}})\f$. Row indices correspond to subjects and
 * column indices to classes.
 */
const double * SubjectClassDistributions(const ClassAssignment *);

/**
 * Returns the estimated number of subjects in each class given this assignment.
 * In Winsberg and De Soete 1993, this vector is represented by \f$\mathbf{F}\f$.
 */
const double * ClassSizes(const ClassAssignment *);

/**
 * Returns the marginalised distribution function over all classes. In Winsberg
 * and De Soete 1993, this vector is represented by 
 * \f$\boldsymbol{\hat{\lambda}}\f$.
 */
const double * UnconditionalClassDistribution(const ClassAssignment *);
