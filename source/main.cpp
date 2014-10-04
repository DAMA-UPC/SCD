/*SCD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <common/types.h>
#include <common/time.h>
#include <communities/communities.h>
#include <graph/graph.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define CHECK_ARGUMENT_STRING(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = argv[index+1]; \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_FLOAT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atof(argv[index+1]); \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_INT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atoi(argv[index+1]); \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_FLAG(index, option,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
		}

using namespace scd;

static void PrintUsage() {
    printf("Usage: scd <flags>\n");
    printf("Availaible flags:\n");
    printf("\t-f [network file name] : Specifies the network file.\n");
    printf("\t-o [output file name] : Specifies the output file name.\n");
    printf("\t-n [number of threads]: Specifies the number of threads to run the algorithm.\n");
    printf("\t-l [lookahead size]: Sets the size of the lookahead iterations to look.\n");
    printf("\t-p [partition file name]: Specifies the partition file name to start the refinement from.\n");
    printf("\t-a [lookahead size]: Specifies the alfa parameter to control the level of cohesion of the communities. Default: 1.0.\n");
}


int main(int argc, char ** argv) {

    bool graphFileNameSet = false;
    bool outputFileNameSet = false;
    bool partitionFileNameSet = false;
    bool numThreadsSet = false;
    bool alfaSet = false;
    char_t * graphFileName = NULL;
    char_t * outputFileName = NULL;
    char_t * partitionFileName = NULL;
    uint32_t numThreads = omp_get_num_procs();
    uint32_t lookahead = 5;
    bool lookaheadSet = false;
    double64_t alfa = 1.0;

    for (uint32_t i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-p", partitionFileName, partitionFileNameSet)
        CHECK_ARGUMENT_INT(i, "-n", numThreads, numThreadsSet)
        CHECK_ARGUMENT_INT(i, "-l", lookahead, lookaheadSet)
        CHECK_ARGUMENT_FLOAT(i, "-a", alfa, alfaSet)
    }

    if (!graphFileNameSet) {
        printf("Graph filename not set\n");
        PrintUsage();
        return 1;
    }
    
    if (numThreads <= 0) {
        printf("Invalid number of threads\n");
        PrintUsage();
        return 2;
    }

    if (!outputFileNameSet) {
        outputFileName = new char_t[512];
        sprintf(outputFileName, "communities.dat");
    }

    CGraph graph;
    uint64_t totalTime = 0,
             initTime = 0, 
             spentTime = 0, 
             loadingTime = 0,
             algorithmTime = 0;


    //==================== LOAD THE GRAPH ==================================
    initTime = StartClock();
    printf("Graph: %s\n", graphFileName);
    printf("OutputFile: %s\n", outputFileName);
    graph.Load(graphFileName, numThreads);
    spentTime = StopClock(initTime);
    loadingTime = spentTime;
    totalTime += spentTime;
    printf("Load time: %lu ms\n", spentTime);
    //======================================================================


    //================ REMOVE EDGES WITHOUT TRIANGLES ======================
    initTime = StartClock();
    printf("Removing edges without triangles ...\n");
    graph.RemoveEdgesNoTriangles(numThreads);
    spentTime = StopClock(initTime);
    algorithmTime += spentTime;
    totalTime += spentTime;
    printf("Removing edges without triangles time: %lu ms\n", spentTime);
    //======================================================================


    //=================== INITIALIZE PARTITION ============================
    initTime = StartClock();
    CommunityPartition partition;
    if(partitionFileNameSet) {
        printf("Loading partition file %s ... \n", partitionFileName);
        if( LoadPartition(&graph,&partition,partitionFileName, alfa) ) {
            printf("Error loading partition\n");
            return 1;
        }
    } else {
        printf("Initial partition file not set. Computing initial partition ...\n");
        if (InitializeSimplePartition(&graph, &partition, alfa)) {
            printf("Error computing initial partition\n");
            return 1;
        }
    }
    spentTime  = StopClock(initTime);
    totalTime += spentTime;
    printf("Initial partition time: %lu ms\n", spentTime);
    //======================================================================


    //================ TRANSFER NODES AMONG PARTITIONS =====================
    initTime = StartClock();
    if (ImproveCommunities(&graph, &partition, numThreads, lookahead, alfa)) {
        printf("Error while improving communities\n");
        return 1;
    }
    spentTime = StopClock(initTime);
    algorithmTime += spentTime;
    totalTime += spentTime;
    printf("Improvement execution time: %lu ms\n", spentTime);
    //======================================================================


    //======================== PRINT RESULTS ===============================
    initTime = StartClock();
    PrintPartition(&graph, &partition, outputFileName);
    spentTime = StopClock(initTime);
    totalTime += spentTime;
    printf("Print partition time: %lu ms\n", spentTime);
    //======================================================================


    printf("\n");
    printf("\n");
    printf("*******************************************************\n");
    printf("%-32s %-10d\n", "Number of Communities:", partition.m_NumCommunities);
    printf("%-32s %-10f\n", "WCC:", partition.m_WCC / (float32_t) graph.GetNumNodes());
    printf("%-32s %-10lu ms\n", "Loading time:", loadingTime);
    printf("%-32s %-10lu ms\n", "Algorithm time:", algorithmTime);
    printf("%-32s %-10lu ms\n", "Total execution time:", totalTime);
    printf("*******************************************************\n");

    FreeResources(&partition);

    if (!outputFileNameSet) {
        delete [] outputFileName;
    }
    return 0;
}



