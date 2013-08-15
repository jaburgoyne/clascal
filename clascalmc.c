// Monte-Carlo tests for CLASCAL
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
#include <time.h>
#include <unistd.h>
#include "inlines.h"
#include "StimulusSet.h"
#include "SubjectSet.h"
#include "Model.h"
#include "Experiment.h"
#include "ClassAssignment.h"
#include "ModelSpace.h"
#include "Solution.h"

static const size_t DEFAULT_SAMPLE_SIZE = 100;

int main(int argc, char * argv[])
{
        // The following calls are necessary for multi-threaded LAPACK. 
        // (void) slamch_("e");
        // (void) dlamch_("e");
        const char * restrict args;
        args = "g:n:";
        int option;
        size_t sampleSize = DEFAULT_SAMPLE_SIZE;
        srand((unsigned)time(NULL));
        while ((option = getopt(argc, argv, args)) != -1) {
                switch (option) {
                        case 'g':
                                if (!optarg) goto Usage;
                                srand((unsigned)strtoul(optarg, NULL, 10));
                                break;
                        case 'n':
                                if (!optarg) goto Usage;
                                sampleSize = SizeRead(optarg);
                                break;
                        default:
                                goto Usage;
                }
        }
        if (argc - optind != 3) goto Usage;
        char * experimentFilename = argv[optind++];
        if (access(experimentFilename, R_OK) != 0) goto Usage;
        Experiment * experiment;
        experiment = NewExperimentFromFilename(experimentFilename);
        char * nullFilename = argv[optind++];
        if (access(nullFilename, R_OK) != 0) goto Usage;
        char * altFilename = argv[optind];
        if (access(altFilename, R_OK) != 0) goto Usage;
        Solution * nullSolution;
        nullSolution = NewSolutionFromFilename(nullFilename, experiment);
        if (!nullSolution)
                ExitWithError("Could not read null hypothesis file");
        Solution * altSolution;
        altSolution = NewSolutionFromFilename(altFilename, experiment);
        if (!altSolution)
                ExitWithError("Could not read alternative hypothesis file");
        const ModelSpace * restrict nullSpace;
        nullSpace = SolutionModelSpace(nullSolution);
        const Model * restrict nullModel = ModelSpaceModel(nullSpace);
        const ModelSpace * restrict altSpace = SolutionModelSpace(altSolution);
        const Model * restrict altModel = ModelSpaceModel(altSpace);
        if (DimensionCount(nullModel) > DimensionCount(altModel)
            || ClassCount(nullModel) > ClassCount(altModel)
            || (ModelSpecificityType(nullModel) == ClassSpecificities
                && ModelSpecificityType(altModel) != ClassSpecificities)
            || (ModelSpecificityType(nullModel) == GlobalSpecificities                
                && ModelSpecificityType(altModel) == NoSpecificities)
            || (DimensionCount(nullModel) == 0
                && ModelSpecificityType(nullModel) == NoSpecificities
                && (DimensionCount(altModel) > 0 
                    || ModelSpecificityType(altModel) != NoSpecificities)))
                ExitWithError("Models are not nested");
        const double logLikelihoodRatio = (2.0 
                                           * (LogLikelihood(altSolution)
                                              - LogLikelihood(nullSolution)));
        const Parameters * nullParameters = SolutionParameters(nullSolution);
        const Parameters * altParameters = SolutionParameters(altSolution);
        size_t lesserRatioCount = 0;
        size_t greaterRatioCount = 0;
        for (size_t n = 0; n < sampleSize; n++) {
                Experiment * restrict sampleExperiment;
                sampleExperiment = NewMonteCarloExperiment(nullSolution);
                Solution * restrict nullFit;
                // Winsberg forces k-means or INDSCAL starts and random initial
                // initial coordinates in most cases.
                nullFit = NewSolutionForExperimentAndModel(sampleExperiment, 
                                                           nullModel, 
                                                           nullParameters);
                Solution * restrict altFit;
                altFit = NewSolutionForExperimentAndModel(sampleExperiment, 
                                                          altModel, 
                                                          altParameters);
                const double sampleLogRatio = (2.0 
                                               * (LogLikelihood(altFit)
                                                  - LogLikelihood(nullFit)));
                if (isless(sampleLogRatio, logLikelihoodRatio)) 
                        lesserRatioCount++;
                if (isgreater(sampleLogRatio, logLikelihoodRatio)) 
                        greaterRatioCount++;
                DeleteSolution(altFit);
                DeleteSolution(nullFit);
                DeleteExperiment(sampleExperiment);
        }
        printf("{\n\t\"lowerBound\": %f,\n\t\"upperBound\": %f\n}\n",
               (double)greaterRatioCount / (double)SizeSum(sampleSize, 1),
               (1.0 
                - (double)lesserRatioCount / (double)SizeSum(sampleSize, 1)));
        DeleteSolution(altSolution);
        DeleteSolution(nullSolution);
        DeleteExperiment(experiment);
        exit(EXIT_SUCCESS);
Usage:
        fprintf(stderr, 
                "Usage: clascalmc [-g uint] [−n uint] exp sol0 sol1\n");
        exit(EXIT_FAILURE);
}
