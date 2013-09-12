// CLASCAL (ExperimentFile.c)
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
#include "inlines.h"
#include "jansson.h"
#include "StimulusSet.h"
#include "_StimulusSet.h"
#include "SubjectSet.h"
#include "Model.h"
#include "Experiment.h"

static Experiment * 
NewExperimentFromJSON(const json_t * restrict experimentJSON)
{
        // There are few checks for type in this function because Jansson
        // returns NULL on error and NULL parameters simply yield a NULL 
        // Experiment.
        const json_t * restrict descriptionJSON = NULL;
        descriptionJSON = json_object_get(experimentJSON, "Description");
        const char * restrict description = json_string_value(descriptionJSON);
        char * restrict descriptionCopy = NULL;
        if (description) {
                size_t descriptionLength = SizeSum(strlen(description), 1);
                descriptionCopy = SafeMalloc(descriptionLength, sizeof(char));
                strcpy(descriptionCopy, description);
        }
        const json_t * restrict stimuliJSON = NULL;
        stimuliJSON = json_object_get(experimentJSON, "StimulusNames");
        const size_t stimulusCount = json_array_size(stimuliJSON);
        if (stimulusCount < 2)
                ExitWithError("Experiment lacks stimulus names");
        char * * restrict stimulusNames = NULL;
        stimulusNames = SafeMalloc(stimulusCount, sizeof(char *));
        for (size_t j = 0; j < stimulusCount; j++) {
                const json_t * restrict stimulusNameJSON = NULL;
                stimulusNameJSON = json_array_get(stimuliJSON, (unsigned int)j);
                const char * restrict stimulusName;
                stimulusName = json_string_value(stimulusNameJSON);
                if (!stimulusName) 
                        ExitWithError("Stimulus names must be strings");
                size_t nameLength = SizeSum(strlen(stimulusName), 1);
                stimulusNames[j] = SafeMalloc(nameLength, sizeof(char));
                strcpy(stimulusNames[j], stimulusName);
        }
        StimulusSet * restrict stims = NULL;
        stims = NewStimulusSet(stimulusCount, stimulusNames);
        json_t * restrict subjectsJSON = NULL;
        subjectsJSON = json_object_get(experimentJSON, "SubjectData");
        const size_t subjectCount = json_object_size(subjectsJSON);
        if (!subjectCount)
                ExitWithError("Experiment lacks named dissimilarity data");
        char * * restrict subjectNames = NULL;
        subjectNames = SafeMalloc(subjectCount, sizeof(char *));
        double * restrict dissimilarities = NULL;
        const size_t pCount = StimulusPairCount(stims);
        const size_t dissimilaritiesSize = SizeProduct(subjectCount, pCount);
        dissimilarities = SafeMalloc(dissimilaritiesSize, sizeof(double));
        void * iterator = json_object_iter(subjectsJSON);
        for (size_t i = 0; i < subjectCount; i++) {
                const char * restrict subjectName = NULL;
                subjectName = json_object_iter_key(iterator);
                if (!subjectName)
                        ExitWithError("Subject names must be strings");
                size_t nameLength = SizeSum(strlen(subjectName), 1);
                subjectNames[i] = SafeMalloc(nameLength, sizeof(char));
                strcpy(subjectNames[i], subjectName);
                const json_t * restrict subjectJSON = NULL;
                subjectJSON = json_object_iter_value(iterator);
                const size_t lineCount = json_array_size(subjectJSON);
                const json_t * restrict firstItem = NULL;
                firstItem = json_array_get(subjectJSON, 0);
                if (json_is_array(firstItem)) {
                        // Skip the first line if given the whole array and
                        // ignore diagonal elements.
                        const size_t off = (lineCount >= stimulusCount) ? 1 : 0;
                        size_t m = 0;
                        for (size_t j = 0; j < stimulusCount - 1; j++) {
                                const json_t * restrict line = NULL;
                                line = json_array_get(subjectJSON,
                                                      (unsigned int)(j + off));
                                for (size_t k = 0; k < j + 1; k++) {
                                        const json_t * restrict datum = NULL;
                                        datum = json_array_get(line,
                                                               (unsigned int)k);
                                        double d = NAN;
                                        if (json_is_real(datum)) {
                                                d = json_real_value(datum);
                                        } else if (json_is_integer(datum)) {
                                                d = (double)json_integer_value(datum);
                                        } else if (json_is_null(datum)) {
                                                d = NAN;
                                        } else {
                                                ExitWithError("Some"
                                                              " dissimilarities"
                                                              " are invalid");
                                        }
                                        dissimilarities[pCount * i + m++] = d;
                                }
                        }
                } else if (json_is_object(firstItem)) {
                        for (size_t m = 0; m < pCount; m++)
                                dissimilarities[pCount * i + m] = NAN;
                        for (size_t l = 0; l < lineCount; l++) {
                                const json_t * restrict datum = NULL;
                                datum = json_array_get(subjectJSON,
                                                       (unsigned int)l);
                                StimulusPair p = {SIZE_MAX, SIZE_MAX};
                                const json_t * restrict j = NULL;
                                j = json_object_get(datum, "j");
                                const json_t * restrict k = NULL;
                                k = json_object_get(datum, "k");
                                if json_is_integer(j)
                                        p.j = (size_t)json_integer_value(j);
                                else ExitWithError("Invalid j-coordinate");
                                if json_is_integer(k)
                                        p.k = (size_t)json_integer_value(k);
                                else ExitWithError("Invalid k-coordinate");
                                size_t m = PairNumberForStimulusPair(stims, &p);
                                const json_t * restrict v = NULL;
                                v = json_object_get(datum, "d");
                                double d = NAN;
                                if (json_is_real(v)) {
                                        d = json_real_value(v);
                                } else if (json_is_integer(v)) {
                                        d = (double)json_integer_value(v);
                                } else if (json_is_null(v)) {
                                        d = NAN;
                                } else {
                                        ExitWithError("Some"
                                                      " dissimilarities"
                                                      " are invalid");
                                }
                                dissimilarities[pCount * i + m] = d;
                        }
                }
                iterator = json_object_iter_next(subjectsJSON, iterator);
        }
        SubjectSet * restrict subjects = NULL;
        subjects = NewSubjectSet(subjectCount, subjectNames);
        return NewExperiment(descriptionCopy, 
                             stims, 
                             subjects,
                             dissimilarities);
}

Experiment * NewExperimentFromFile(FILE * restrict f)
{
        json_error_t jsonError;
        json_t * restrict experimentJSON = json_loadf(f, &jsonError);
        if (!experimentJSON) {
                fprintf(stderr,
                        "WARNING: JSON parsing failed in line %i (%s).\n",
                        jsonError.line,
                        jsonError.text);
                return NULL;
        }
        Experiment * experiment = NewExperimentFromJSON(experimentJSON);
        json_decref(experimentJSON);
        return experiment;
}

Experiment * NewExperimentFromFilename(const char * restrict filename)
{
        json_error_t jsonError;
        json_t * restrict experimentJSON = json_load_file(filename, &jsonError);
        if (!experimentJSON) {
                fprintf(stderr,
                        "WARNING: JSON parsing failed in line %i of %s (%s).\n",
                        jsonError.line,
                        filename,
                        jsonError.text);
                return NULL;
        }
        Experiment * experiment = NewExperimentFromJSON(experimentJSON);
        json_decref(experimentJSON);
        return experiment;
}
