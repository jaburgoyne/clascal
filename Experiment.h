// CLASCAL (Experiment.h)
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
 * @file Experiment.h
 *
 * Experiment type and related definitions. Requires SubjectSet.h and 
 * StimulusSet.h.
 */

struct _Experiment;
/**
 * An experiment involving dissimilarity ratings.
 */
typedef struct _Experiment Experiment;

/**
 * Allocates, initialises, and returns a new Experiment.
 *
 * @param description  a description of the experiment (may be NULL)
 * @param stimulusSet  the stimulus set used for the experiment
 * @param subjectSet  the subject set used for the experiment
 * @param dissimilarities  dissimilarity data for the experiment, row indices
 *                         corresponding to subjects and column indices to
 *                         stimulus pairs
 */
Experiment * NewExperiment(char *,
                           StimulusSet *,
                           SubjectSet *,
                           double *);

/**
 * Allocates, initialises, and returns a new Experiment from a file.
 *
 * @param file  the file where the Experiment data are located
 */
Experiment * NewExperimentFromFile(FILE *);

/**
 * Allocates, initialises, and returns a new Experiment from a filename.
 *
 * @param filename  the path of the file where the Experiment data located
 */
Experiment * NewExperimentFromFilename(const char *);

/**
 * Deallocates an Experiment.
 */
void DeleteExperiment(Experiment *);

/**
 * Returns a description of this experiment.
 */
const char * ExperimentDescription(const Experiment *);

/**
 * Returns the StimulusSet used for an experiment.
 */
const StimulusSet * ExperimentStimulusSet(const Experiment *);

/**
 * Returns the SubjectSet used for an experiment.
 */
const SubjectSet * ExperimentSubjectSet(const Experiment *);

/**
 * Returns the dissimilarity ratings from the experiment. In Winsberg and De
 * Soete, this matrix is represented by \f$Y\f$. Row indices correspond
 * to subjects and column indices to stimulus pair numbers.
 */
const double * SubjectDissimilarities(const Experiment *);

