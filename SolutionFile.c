// CLASCAL (SolutionFile.c)
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
#include "SubjectSet.h"
#include "Experiment.h"
#include "Model.h"
#include "ClassAssignment.h"
#include "ModelSpace.h"
#include "_ModelSpace.h"
#include "Solution.h"
#include "_Solution.h"

static json_t * json_array_with_vector(const double * restrict vector, 
                                       size_t length, 
                                       size_t stride)
{
        json_t * restrict array = json_array();
        for (size_t i = 0; i < SizeProduct(length, stride); i += stride)
                json_array_append_new(array, json_real(vector[i]));
        return array;
}

/* N.B. This function leaks memory when loading the parameters and the model. */
static Solution * NewSolutionFromJSON(const json_t * restrict solutionJSON,
                                      const Experiment * restrict experiment)
{
        // There are few checks for type in this function because Jansson
        // returns NULL on error and NULL parameters simply yield a NULL
        // Solution.
        const StimulusSet * restrict stimulusSet = NULL;
        stimulusSet = ExperimentStimulusSet(experiment);
        const json_t * restrict coordinatesJSON = NULL;
        coordinatesJSON = json_object_get(solutionJSON, "Coordinates");
        const json_t * restrict stimulusJSON = NULL;
        stimulusJSON = json_object_get(coordinatesJSON,
                                       StimulusNameForIndex(stimulusSet, 0));
        const size_t dimensionCount = json_array_size(stimulusJSON);
        const SubjectSet * restrict subjectSet = NULL;
        subjectSet = ExperimentSubjectSet(experiment);
        const json_t * restrict posteriorsJSON = NULL;
        posteriorsJSON = json_object_get(solutionJSON, 
                                         "PosteriorProbabilities");
        const json_t * restrict subjectJSON = NULL;
        subjectJSON = json_object_get(posteriorsJSON, 
                                      SubjectNameForIndex(subjectSet, 0));
        const size_t classCount = json_array_size(subjectJSON);
        const json_t * restrict specificitiesJSON = NULL;
        specificitiesJSON = json_object_get(solutionJSON, "Specificities");
        const json_t * restrict specJSON = NULL;
        specJSON = json_object_get(specificitiesJSON, 
                                   StimulusNameForIndex(stimulusSet, 0));
        const SpecificityType specType = (specificitiesJSON
                                          ? (json_is_number(specJSON)
                                             ? GlobalSpecificities
                                             : ClassSpecificities)
                                          : NoSpecificities);
        // N.B. This Model will never be deleted.
        Model * restrict model = NewModel(dimensionCount, classCount, specType);
        const size_t subjectCount = SubjectCount(subjectSet);
        const size_t distributionsSize = SizeProduct(subjectCount, classCount);
        double * restrict distributions = NULL;
        distributions = SafeMalloc(distributionsSize, sizeof(double));
        for (size_t i = 0; i < subjectCount; i++) {
                const char * subjectName = SubjectNameForIndex(subjectSet, i);
                subjectJSON = json_object_get(posteriorsJSON, subjectName);
                for (size_t t = 0; t < classCount; t++) {
                        const json_t * restrict valueJSON = NULL;
                        valueJSON = json_array_get(subjectJSON, t);
                        const double value = json_number_value(valueJSON);
                        distributions[classCount * i + t] = value;
                }
        }
        ClassAssignment * restrict assignment = NULL;
        assignment = NewClassAssignment(subjectSet, model, distributions);
        free(distributions);
        const size_t weightsLength = ((specType == GlobalSpecificities)
                                      ? SizeSum(dimensionCount, 1)
                                      : dimensionCount);
        const size_t weightsSize = SizeProduct(weightsLength, classCount);
        double * restrict weights = NULL;
        if (classCount > 1 && (dimensionCount > 0 || specType)) {
                weights = SafeMalloc(weightsSize, sizeof(double));
                const json_t * restrict weightsJSON = NULL;
                weightsJSON = json_object_get(solutionJSON, "Weights");
                for (size_t t = 0; t < classCount; t++) {
                        const json_t * restrict classJSON = NULL;
                        classJSON = json_array_get(weightsJSON, t);
                        for (size_t r = 0; r < weightsLength; r++) {
                                const json_t * restrict valJSON = NULL;
                                valJSON = json_array_get(classJSON, r);
                                const double value = json_number_value(valJSON);
                                weights[classCount * r + t] = value;
                        }
                }
        }
        const size_t stimulusCount = StimulusCount(stimulusSet);
        const size_t coordinatesSize = SizeProduct(stimulusCount, 
                                                   dimensionCount); 
        double * restrict coordinates = NULL;
        if (dimensionCount > 0) {
                coordinates = SafeMalloc(coordinatesSize, sizeof(double));
                for (size_t j = 0; j < stimulusCount; j++) {
                        const char * stName = NULL;
                        stName = StimulusNameForIndex(stimulusSet, j);
                        stimulusJSON = json_object_get(coordinatesJSON, stName);
                        for (size_t r = 0; r < dimensionCount; r++) {
                                const json_t * restrict valJSON = NULL;
                                valJSON = json_array_get(stimulusJSON, r);
                                const double value = json_number_value(valJSON);
                                coordinates[dimensionCount * j + r] = value;
                        }
                }
        }
        double * restrict specificities = NULL;
        switch (specType) {
                case ClassSpecificities:
                        specificities = SafeMalloc(SizeProduct(classCount,
                                                               stimulusCount), 
                                                   sizeof(double));
                        for (size_t j = 0; j < stimulusCount; j++) {
                                const char * stimName = NULL;
                                stimName = StimulusNameForIndex(stimulusSet, j);
                                specJSON = json_object_get(specificitiesJSON, 
                                                           stimName);
                                for (size_t t = 0; t < classCount; t++) {
                                        const json_t * restrict valJSON = NULL;
                                        valJSON = json_array_get(specJSON, t);
                                        const double value 
                                        = json_number_value(valJSON);
                                        specificities[stimulusCount * t + j] 
                                        = value;
                                }
                        }
                        break;
                case GlobalSpecificities:
                        specificities = SafeMalloc(stimulusCount, 
                                                   sizeof(double));
                        for (size_t j = 0; j < stimulusCount; j++) {
                                const char * stimName = NULL;
                                stimName = StimulusNameForIndex(stimulusSet, j);
                                specJSON = json_object_get(specificitiesJSON, 
                                                           stimName);
                                const double value 
                                = json_number_value(specJSON);
                                specificities[j] = value;
                        }
                        break;
                default:
                        // Leave specificities NULL.
                        break;
        }
        const size_t pairCount = StimulusPairCount(stimulusSet);
        const size_t distancesSize = SizeProduct(subjectCount, pairCount);
        ModelSpace * space = NULL;
        if (dimensionCount || specificities) {
                space = NewModelSpace(stimulusSet, 
                                      model, 
                                      weights, 
                                      coordinates, 
                                      specificities);
        } else {
                double * distances = SafeMalloc(distancesSize, sizeof(double));
                const json_t * restrict dissJSON = NULL;
                dissJSON = json_object_get(solutionJSON, 
                                           "PredictedDissimilarities");
                for (size_t t = 0; t < classCount; t++) {
                        const json_t * restrict classJSON = NULL;
                        classJSON = json_array_get(dissJSON,
                                                   (unsigned int)t);
                        size_t m = 0;
                        for (size_t j = 0; j < stimulusCount - 1; j++) {
                                const json_t * restrict line = NULL;
                                line = json_array_get(classJSON,
                                                      (unsigned int)j);
                                for (size_t k = 0; k < j + 1; k++) {
                                        const json_t * restrict datum = NULL;
                                        datum = json_array_get(line, 
                                                               (unsigned int)k);
                                        if (!json_is_number(datum))
                                                ExitWithError("Some"
                                                              " dissimilarities"
                                                              " are missing");
                                        double d = json_real_value(datum);
                                        distances[pairCount * t + m++] = d;
                                }
                        }
                }
                space = NewDummySpace(stimulusSet, model, distances);
        }
        free(specificities);
        free(coordinates);
        free(weights);
        const json_t * restrict priorJSON = NULL;
        priorJSON = json_object_get(solutionJSON, "PriorDistribution");
        double * restrict prior = SafeMalloc(classCount, sizeof(double));
        for (size_t t = 0; t < classCount; t++)
                prior[t] = json_number_value(json_array_get(priorJSON, t));
        const json_t * restrict estimatedVarianceJSON = NULL;
        estimatedVarianceJSON = json_object_get(solutionJSON,
                                                "EstimatedVariance");
        const double sigmaSquared = json_number_value(estimatedVarianceJSON);
        const json_t * restrict parametersJSON = NULL;
        parametersJSON = json_object_get(solutionJSON, "Parameters");
        // N.B. This Parameters object will never be freed.
        Parameters * parameters = SafeMalloc(1, sizeof(Parameters));
        // Verbosity is not saved in the file and defaults to zero.
        parameters->verbosity = SILENT;
        const json_t * restrict useTorgersonStartJSON = NULL;
        useTorgersonStartJSON = json_object_get(parametersJSON,
                                                "TorgersonStart");
        parameters->useTorgersonStart = (json_is_true(useTorgersonStartJSON)
                                         ? true
                                         : false);
        const json_t * restrict useKMeansStartJSON = NULL;
        useKMeansStartJSON = json_object_get(parametersJSON,
                                             "KMeansStart");
        parameters->useKMeansStart = (json_is_true(useKMeansStartJSON)
                                      ? true
                                      : false);
        parameters->maxEMIterationCount 
        = json_integer_value(json_object_get(parametersJSON,
                                             "MaxEMIterationCount"));
        parameters->EMImprovementAmount 
        = json_real_value(json_object_get(parametersJSON,
                                          "EMImprovementAmount"));
        parameters->maxMIterationCount 
        = json_integer_value(json_object_get(parametersJSON,
                                             "MaxMIterationCount"));
        parameters->MImprovementFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "MImprovementFactor"));
        parameters->maxSpatialIterationCount 
        = json_integer_value(json_object_get(parametersJSON,
                                             "MaxSpatialIterationCount"));
        parameters->minSpecificityValue 
        = json_real_value(json_object_get(parametersJSON,
                                          "MinSpecificityValue"));        
        parameters->spatialImprovementFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "SpatialImprovementFactor"));
        parameters->maxWeightIterationCount 
        = json_integer_value(json_object_get(parametersJSON,
                                             "MaxWeightIterationCount"));
        parameters->minWeightValue 
        = json_real_value(json_object_get(parametersJSON,
                                          "MinWeightValue"));        
        parameters->weightImprovementFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "WeightImprovementFactor"));
        parameters->maxSearchIterationCount 
        = json_integer_value(json_object_get(parametersJSON,
                                             "MaxSearchIterationCount"));
        parameters->firstStepSize 
        = json_real_value(json_object_get(parametersJSON,
                                          "FirstStepSize"));
        parameters->margin 
        = json_real_value(json_object_get(parametersJSON,
                                          "Margin"));
        parameters->invConditionNumber 
        = json_real_value(json_object_get(parametersJSON,
                                          "InverseConditionNumber"));
        parameters->slopeImprovementFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "SlopeImprovementFactor"));
        parameters->stepIncreaseFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "StepIncreaseFactor"));
        parameters->stepDecreaseFactor 
        = json_real_value(json_object_get(parametersJSON,
                                          "StepDecreaseFactor"));
        Solution * solution = NewSolution(experiment,
                                          assignment,
                                          space,
                                          prior,
                                          sigmaSquared,
                                          parameters);
        free(prior);
        return solution;
}

Solution * NewSolutionFromFile(FILE * restrict f, 
                               const Experiment * restrict experiment)
{
        json_error_t jsonError;
        json_t * restrict solutionJSON = json_loadf(f, &jsonError);
        if (!solutionJSON) return NULL;
        Solution * solution = NewSolutionFromJSON(solutionJSON, experiment);
        json_decref(solutionJSON);
        return solution;
}

Solution * NewSolutionFromFilename(const char * restrict filename,
                                   const Experiment * restrict experiment)
{
        json_error_t jsonError;
        json_t * restrict solutionJSON = json_load_file(filename, &jsonError);
        if (!solutionJSON) return NULL;
        Solution * solution = NewSolutionFromJSON(solutionJSON, experiment);
        json_decref(solutionJSON);
        return solution;
}

static json_t * SolutionJSON(const Solution * restrict self)
{
        json_t * restrict root = json_object();
        const ClassAssignment * restrict assignment = NULL;
        assignment = SolutionClassAssignment(self);
        const Model * restrict model = ClassAssignmentModel(assignment);
        const size_t classCount = ClassCount(model);
        const SubjectSet * restrict subjectSet = NULL;
        subjectSet = ClassAssignmentSubjectSet(assignment);
        const size_t subjectCount = SubjectCount(subjectSet);
        const double * restrict subjectClassDistributions = NULL;
        subjectClassDistributions = SubjectClassDistributions(assignment);
        json_t * restrict subjectClassDistributionsJSON = json_object();
        for (size_t i = 0; i < subjectCount; i++) {
                const double * distribution = NULL;
                distribution = subjectClassDistributions + classCount * i;
                json_t * restrict distributionJSON = NULL;
                distributionJSON = json_array_with_vector(distribution, 
                                                          classCount,
                                                          1);
                const char * name = SubjectNameForIndex(subjectSet, i);
                json_object_set_new(subjectClassDistributionsJSON, 
                                    name ? name : "", 
                                    distributionJSON);
        }
        json_object_set_new(root, 
                            "PosteriorProbabilities",
                            subjectClassDistributionsJSON);
        const ModelSpace * restrict space = SolutionModelSpace(self);
        const size_t dimensionCount = DimensionCount(model);
        const SpecificityType specificityType = ModelSpecificityType(model);
        const StimulusSet * restrict stimulusSet = NULL;
        stimulusSet = ModelSpaceStimulusSet(space);
        const size_t stimulusCount = StimulusCount(stimulusSet);
        if (classCount > 1 && (dimensionCount > 0 || specificityType)) {
                const size_t weightsLength = ((specificityType 
                                               == GlobalSpecificities)
                                              ? SizeSum(dimensionCount, 1)
                                              : dimensionCount);
                const double * restrict weights = Weights(space);
                json_t * restrict weightsJSON = json_array();
                for (size_t t = 0; t < classCount; t++) {
                        const double * classWeights = weights + t;
                        json_t * restrict classWeightsJSON = NULL;
                        classWeightsJSON = json_array_with_vector(classWeights, 
                                                                  weightsLength, 
                                                                  classCount);
                        json_array_append_new(weightsJSON, classWeightsJSON);
                }
                json_object_set_new(root, "Weights", weightsJSON);
        }
        if (dimensionCount > 0) {
                const double * restrict coordinates = Coordinates(space);
                json_t * restrict coordinatesJSON = json_object();
                for (size_t j = 0; j < stimulusCount; j++) {
                        const double * row = coordinates + dimensionCount * j;
                        json_t * restrict rowJSON = NULL;
                        rowJSON = json_array_with_vector(row, 
                                                         dimensionCount, 
                                                         1);
                        const char * name = NULL;
                        name = StimulusNameForIndex(stimulusSet, j);
                        json_object_set_new(coordinatesJSON, 
                                            name ? name : "", rowJSON);
                }
                json_object_set_new(root, "Coordinates", coordinatesJSON);
        }
        if (specificityType) {
                const double * specificities = Specificities(space);
                json_t * restrict specificitiesJSON = json_object();
                for (size_t j = 0; j < stimulusCount; j++) {
                        json_t * restrict spJSON = NULL;
                        if (specificityType == ClassSpecificities) {
                                const double * col = specificities + j;
                                spJSON = json_array_with_vector(col, 
                                                                classCount, 
                                                                stimulusCount);
                        } else {
                                spJSON = json_real(specificities[j]);
                        }
                        const char * name = NULL;
                        name = StimulusNameForIndex(stimulusSet, j);
                        json_object_set_new(specificitiesJSON, 
                                            name ? name : "", 
                                            spJSON);
                }
                json_object_set_new(root, 
                                    "Specificities", 
                                    specificitiesJSON);
        }
        json_object_set_new(root, 
                            "ParameterCount",
                            json_integer((int)ParameterCount(space)));
        json_object_set_new(root,
                            "DegreesOfFreedom",
                            json_integer((int)DegreesOfFreedom(self)));
        json_t * restrict distancesJSON = json_array();
        for (size_t t = 0; t < classCount; t++) {
                json_t * restrict classMeansJSON = json_array();
                for (size_t k = 1; k < stimulusCount; k++) {
                        json_t * restrict rowJSON = json_array();
                        for (size_t j = 0; j < k; j++) {
                                const StimulusPair pair = { j, k };
                                const double d = DistanceForPairAndClass(space, 
                                                                         &pair,
                                                                         t);
                                json_array_append_new(rowJSON, json_real(d));
                        }
                        json_array_append_new(classMeansJSON, rowJSON);
                }
                json_array_append_new(distancesJSON, classMeansJSON);
        }
        json_object_set_new(root, "PredictedDissimilarities", distancesJSON);
        const double * prior = PriorDistribution(self);
        json_t * restrict priorJSON = json_array();
        for (size_t t = 0; t < classCount; t++)
                json_array_append_new(priorJSON, json_real(prior[t]));
        json_object_set_new(root, 
                            "PriorDistribution", 
                            priorJSON);
        json_object_set_new(root, 
                            "EstimatedVariance", 
                            json_real(EstimatedVariance(self)));
        const Parameters * restrict parameters = SolutionParameters(self);
        json_t * restrict parametersJSON = json_object();
        json_object_set_new(parametersJSON, 
                            "TorgersonStart", 
                            (parameters->useTorgersonStart 
                             ? json_true() 
                             : json_false()));
        json_object_set_new(parametersJSON, 
                            "KMeansStart", 
                            (parameters->useKMeansStart 
                             ? json_true() 
                             : json_false()));
        json_object_set_new(parametersJSON, 
                            "MaxEMIterationCount", 
                            json_integer((int)parameters->maxEMIterationCount));
        json_object_set_new(parametersJSON, 
                            "EMImprovementAmount",
                            json_real(parameters->EMImprovementAmount));
        json_object_set_new(parametersJSON, 
                            "MaxMIterationCount",
                            json_integer((int)parameters->maxMIterationCount));
        json_object_set_new(parametersJSON, 
                            "MImprovementFactor",
                            json_real(parameters->MImprovementFactor));
        json_object_set_new(parametersJSON, 
                            "MaxSpatialIterationCount",
                            json_integer((int)parameters->maxSpatialIterationCount));
        json_object_set_new(parametersJSON, 
                            "MinSpecificityValue",
                            json_real(parameters->minSpecificityValue));
        json_object_set_new(parametersJSON, 
                            "SpatialImprovementFactor",
                            json_real(parameters->spatialImprovementFactor));
        json_object_set_new(parametersJSON, 
                            "MaxWeightIterationCount",
                            json_integer((int)parameters->maxWeightIterationCount));
        json_object_set_new(parametersJSON, 
                            "MinWeightValue",
                            json_real(parameters->minWeightValue));
        json_object_set_new(parametersJSON, 
                            "WeightImprovementFactor",
                            json_real(parameters->weightImprovementFactor));
        json_object_set_new(parametersJSON, 
                            "MaxSearchIterationCount",
                            json_integer((int)parameters->maxSearchIterationCount));
        json_object_set_new(parametersJSON, 
                            "FirstStepSize",
                            json_real(parameters->firstStepSize));
        json_object_set_new(parametersJSON, 
                            "Margin", 
                            json_real(parameters->margin));
        json_object_set_new(parametersJSON,
                            "InverseConditionNumber",
                            json_real(parameters->invConditionNumber));
        json_object_set_new(parametersJSON, 
                            "SlopeImprovementFactor",
                            json_real(parameters->slopeImprovementFactor));
        json_object_set_new(parametersJSON, 
                            "StepIncreaseFactor",
                            json_real(parameters->stepIncreaseFactor));
        json_object_set_new(parametersJSON, 
                            "StepDecreaseFactor",
                            json_real(parameters->stepDecreaseFactor));
        json_object_set_new(root, "Parameters", parametersJSON);
        json_object_set_new(root, 
                            "LogLikelihood", 
                            json_real(LogLikelihood(self)));
        json_object_set_new(root, 
                            "AkaikeCriterion", 
                            json_real(AkaikeCriterion(self)));
        json_object_set_new(root, 
                            "BayesianCriterion",
                            json_real(BayesianCriterion(self)));
        return root;
}
        
void SaveSolutionToFile(const Solution * restrict self, FILE * restrict f)
{
        json_t * restrict solutionJSON = SolutionJSON(self);
        if (json_dumpf(solutionJSON, f, JSON_INDENT(8)|JSON_SORT_KEYS) != 0)
                ExitWithError("Could not save file");
        json_decref(solutionJSON);
}
            
void SaveSolutionToFilename(const Solution * restrict self,
                            const char * restrict filename)
{
        json_t * restrict solutionJSON = SolutionJSON(self);
        if (0 
            != json_dump_file(solutionJSON, 
                              filename,
                              JSON_INDENT(8)|JSON_SORT_KEYS))
                ExitWithError("Could not save file");
        json_decref(solutionJSON);
}
