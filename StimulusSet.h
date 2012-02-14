// CLASCAL (StimulusSet.h)
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
 * @file StimulusSet.h
 *
 * Stimulus set type and related definitions.
 */

struct _StimulusSet;
/**
 * A set of named stimuli used for an experiment.
 */
typedef struct _StimulusSet StimulusSet;

/**
 * A pair of stimuli.
 */
struct _StimulusPair {
        /** the first member of the pair */
        size_t j;
        /** the second member of the pair */
        size_t k;
};
typedef struct _StimulusPair StimulusPair;

/**
 * Allocates, initialises, and returns a new StimulusSet.
 *
 * @param stimulusCount  the number of stimuli in the set
 * @param stimulusNames  an array of names for each stimulus
 */
StimulusSet * NewStimulusSet(size_t, char * *);

/**
 * Deallocates a StimulusSet.
 */
void DeleteStimulusSet(StimulusSet * restrict self);

/**
 * Returns the number of stimuli in the set. In Winsberg and De Soete 1993, this
 * value is represented by \f$J\f$.
 */
size_t StimulusCount(const StimulusSet *);

/**
 * Returns the number of stimulus pairs in the set. In Winsberg and De Soete
 * 1993, this value is represented by \f$M\f$.
 */
size_t StimulusPairCount(const StimulusSet *);

/**
 * Returns the name of a stimulus given its index within the set.
 */
const char * StimulusNameForIndex(const StimulusSet *, size_t);
