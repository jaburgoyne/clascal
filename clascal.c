// CLASCAL
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

#include <libgen.h>
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

int main(int argc, char * argv[])
{
        const char * restrict args;
        args = "o:r:R:sSt:T:cg:GkC:e:E:j:l:L:m:M:n:N:p:P:v:w:W:x:X:z:";
        int option;
        char workingDirectory[PATH_MAX];
        getcwd(workingDirectory, PATH_MAX);
        char outputDirectory[PATH_MAX];
        getcwd(outputDirectory, PATH_MAX);
        size_t verbosity = 0;
        size_t minR = 0;
        size_t maxR = 0;
        SpecificityType specificityType = NoSpecificities;
        size_t minT = 0;
        size_t maxT = 0;
        srand((unsigned)time(NULL));
        Parameters dereferencedParams = DEFAULT_PARAMETERS;
        Parameters * restrict params = &dereferencedParams;
        while ((option = getopt(argc, argv, args)) != -1) {
                switch (option) {
                        case 'o':
                                if (!optarg 
                                    || strlen(optarg) >= PATH_MAX)
                                        goto Usage;
                                strncpy(outputDirectory, 
                                        optarg, 
                                        PATH_MAX - 1);
                                break;
                        case 'v':
                                if (!optarg) goto Usage;
                                verbosity = SizeRead(optarg);
                                if (verbosity == 0) 
                                        params->verbosity = SILENT;
                                else if (verbosity == 1) 
                                        params->verbosity = VERBOSE;
                                else if (verbosity == 2)
                                        params->verbosity = VERY_VERBOSE;
                                else if (verbosity == 3)
                                        params->verbosity = VERY_VERY_VERBOSE;
                                else params->verbosity = VERY_VERY_VERY_VERBOSE;
                                break;
                        case 'r':
                                if (!optarg) goto Usage;
                                minR = SizeRead(optarg);
                                break;
                        case 'R':
                                if (!optarg) goto Usage;
                                maxR = SizeRead(optarg);
                                break;
                        case 's':
                                // Skip if specificities are already set.
                                if (!specificityType)
                                        specificityType = GlobalSpecificities;
                                break;
                        case 'S':
                                specificityType = ClassSpecificities;
                                break;
                        case 't':
                                if (!optarg) goto Usage;
                                minT = SizeRead(optarg);
                                break;
                        case 'T':
                                if (!optarg) goto Usage;
                                maxT = SizeRead(optarg);
                                break;
                        case 'c':
                                params->useTorgersonStart = true;
                                break;
                        case 'g':
                                if (!optarg) goto Usage;
                                srand((unsigned)strtoul(optarg, NULL, 10));
                                break;
                        case 'G':
                                params->useKMeansStart = false;
                                params->useTorgersonStart = false;
                                break;
                        case 'k':
                                params->useKMeansStart = true;
                                break;
                        case 'C':
                                if (!optarg) goto Usage;
                                params->invConditionNumber 
                                = strtod(optarg, NULL);
                                break;
                        case 'e':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -e set to 1e-3 for"
                                        " compatibility with Winsberg.\n");
#else
                                
                                if (!optarg) goto Usage;
                                params->EMImprovementAmount 
                                = strtod(optarg, NULL);
#endif
                                break;
                        case 'E':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -E set to 150 for"
                                        " compatibility with Winsberg.\n");
#else
                                if (!optarg) goto Usage;
                                params->maxEMIterationCount 
                                = SizeRead(optarg);
#endif
                                break;
                        case 'j':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -j set to 5e-3 for"
                                        " compatibility with Winsberg.\n");
#else                                
                                if (!optarg) goto Usage;
                                params->margin = strtod(optarg, NULL);
#endif
                                break;
                        case 'l':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -l set to 0.1 for"
                                        " compatibility with Winsberg.\n");
#else                                
                                if (!optarg) goto Usage;
                                params->slopeImprovementFactor 
                                = strtod(optarg, NULL);
#endif
                                break;
                        case 'L':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -L set to 5 for"
                                        " compatibility with Winsberg.\n");
#else
                                if (!optarg) goto Usage;
                                params->maxSearchIterationCount
                                = SizeRead(optarg);
#endif
                                break;
                        case 'm':
                                if (!optarg) goto Usage;
                                params->MImprovementFactor
                                = strtod(optarg, NULL);
                                break;
                        case 'M':
                                if (!optarg) goto Usage;
                                params->maxMIterationCount
                                = SizeRead(optarg);
                                break;
                        case 'n':
                                if (!optarg) goto Usage;
                                params->minWeightValue
                                = strtod(optarg, NULL);
                                break;
                        case 'N':
                                if (!optarg) goto Usage;
                                params->minSpecificityValue
                                = strtod(optarg, NULL);
                                break;
                        case 'p':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -p set to 0.5 for"
                                        " compatibility with Winsberg.\n");
#else
                                if (!optarg) goto Usage;
                                params->stepDecreaseFactor
                                = strtod(optarg, NULL);
#endif
                                break;
                        case 'P':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -P set to 2.0 for"
                                        " compatibility with Winsberg.\n");
#else
                                if (!optarg) goto Usage;
                                params->stepIncreaseFactor
                                = strtod(optarg, NULL);
#endif
                                break;
                        case 'w':
                                if (!optarg) goto Usage;
                                params->weightImprovementFactor
                                = strtod(optarg, NULL);
                                break;
                        case 'W':
                                if (!optarg) goto Usage;
                                params->maxWeightIterationCount
                                = SizeRead(optarg);
#ifdef WINSBERG_LEGACY
                                if (params->maxWeightIterationCount > 4) {
                                        params->maxWeightIterationCount = 4;
                                        fprintf(stderr,
                                                "WARNING: Option -W set to 4"
                                                " for compatibility with"
                                                " Winsberg.\n");
                                }
#endif
                                break;
                        case 'x':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -x ignored for"
                                        " compatibility with Winsberg. Will"
                                        " use value of option -w.\n");
#else
                                if (!optarg) goto Usage;
                                params->spatialImprovementFactor
                                = strtod(optarg, NULL);
#endif
                                break;
                        case 'X':
                                if (!optarg) goto Usage;
                                params->maxSpatialIterationCount
                                = SizeRead(optarg);
#ifdef WINSBERG_LEGACY
                                if (params->maxSpatialIterationCount > 7) {
                                        params->maxSpatialIterationCount = 7;
                                        fprintf(stderr,
                                                "WARNING: Option -X set to 7"
                                                " for compatibility with"
                                                " Winsberg.\n");
                                }
#endif
                                break;
                        case 'z':
#ifdef WINSBERG_LEGACY
                                fprintf(stderr,
                                        "WARNING: Option -z set to 1.0 for"
                                        " compatibility with Winsberg.\n");
#else
                                if (!optarg) goto Usage;
                                params->firstStepSize
                                = strtod(optarg, NULL);
#endif
                                break;
                        case '?':
                                fprintf(stderr,
                                        "Unknown or invalid option -%c.\n",
                                        optopt);
                                goto Usage;
                        default:
                                goto Usage;
                }
        }
#ifdef WINSBERG_LEGACY
        params->spatialImprovementFactor = params->weightImprovementFactor;
        if (specificityType == NoSpecificities 
            && params->firstStepSize == 1.0) {
                params->firstStepSize = 0.5;
                fprintf(stderr,
                        "WARNING: Option -z reset to 0.5 for compatibility with"
                        " Winsberg.\n");
        }
#endif
        for ( ; optind < argc; optind++) {
                char * filename = argv[optind];
                if (access(filename, R_OK) != 0) {
                        fprintf(stderr, 
                                "WARNING: '%s' skipped (unreadable).\n", 
                                filename);
                        continue;
                }
                Experiment * experiment = NULL;
                experiment = NewExperimentFromFilename(filename);
                if (maxR < minR) maxR = minR;
                const SubjectSet * subjectSet = NULL;
                subjectSet = ExperimentSubjectSet(experiment); 
                if (minT == 0) minT = SubjectCount(subjectSet);
                if (maxT < minT) maxT = minT;
                for (size_t R = minR; R <= maxR; R++) {
                        for (size_t T = minT; T <= maxT; T++) {
                                if (params->verbosity >= VERBOSE)
                                        fprintf(stdout,
                                                "Fitting to %s"
                                                " with %zu dimensions"
                                                " and %zu classes.\n\n",
                                                filename,
                                                R,
                                                T);
                                Model * model = NewModel(R, T, specificityType);
                                Solution * s = NULL;
                                s = NewSolutionForExperimentAndModel(experiment, 
                                                                     model, 
                                                                     params);
                                char * baseName = basename(filename);
                                char outputName[NAME_MAX];
                                const char * restrict format;
                                switch(specificityType) {
                                        case ClassSpecificities:
                                                format = "T%u-R%u-SS-%s";
                                                break;
                                        case GlobalSpecificities:
                                                format = "T%u-R%u-s-%s";
                                                break;
                                        default:
                                                format = "T%u-R%u-%s";
                                }
                                snprintf(outputName, 
                                         NAME_MAX, 
                                         format,
                                         T,
                                         R,
                                         baseName);
                                chdir(outputDirectory);
                                SaveSolutionToFilename(s, outputName);
                                chdir(workingDirectory);                                
                                DeleteSolution(s);
                                DeleteModel(model);
                        }
                }
                DeleteExperiment(experiment); 
        }
        exit(EXIT_SUCCESS);
Usage:
        fprintf(stderr, 
                "Usage: clascal [−cIksS] [−C real] [−e real] [−E uint]"
                " [−g uint] [-G] [-j real] [−l real] [−L uint] [−m real]"
                " [−M uint] [−n real] [−L uint] [−m real] [−M uint] [−N real]"
                " [−o dir] [−p real] [−P Real] [−r uint] [−R uint] [−t uint]"
                " [−T uint] [−w real] [−W uint] [−x real] [−X uint] [−z real]"
                " file ..."
                "\n");
        exit(EXIT_FAILURE);
}
