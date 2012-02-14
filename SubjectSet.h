// CLASCAL (SubjectSet.h)
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
 * @file SubjectSet.h
 *
 * Subject set type and related definitions.
 */

struct _SubjectSet;
/**
 * A set of subjects (e.g., those used for an experiment).
 */
typedef struct _SubjectSet SubjectSet;

/**
 * Allocated, initialises, and returns a new SubjectSet.
 * 
 * @param subjectCount  the number of subjects in the set
 * @param subjectNames  an array of names (or other identifiers) for each 
 *                      subject
 */
SubjectSet * NewSubjectSet(size_t, char * *);

/**
 * Deallocates a SubjectSet.
 */
void DeleteSubjectSet(SubjectSet * restrict self);

/**
 * Returns the number of subjects in the set. In Winsberg and De Soete 1993, 
 * this value is represented by \f$N\f$.
 */
size_t SubjectCount(const SubjectSet *);

/**
 * Returns the name of a subject given the subject's index.
 */
const char * SubjectNameForIndex(const SubjectSet *, size_t);
