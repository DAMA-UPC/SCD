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
#include <communities/communities.h>
#include <graph/graph.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

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


uint64_t startClock(){
        timeval time;
        gettimeofday(&time, NULL);
        uint64_t initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        return initTime;
}

uint64_t stopClock(uint64_t initTime){
        timeval time;
        gettimeofday(&time, NULL);
        uint64_t endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        return endTime - initTime;
}


int main(int argc, char ** argv) {

	bool graphFileNameSet = false;
	bool outputFileNameSet = false;
	bool numThreadsSet = false;
	char_t * graphFileName = NULL;
	char_t * outputFileName = NULL;
	uint32_t numThreads = 1;

	for( uint32_t i = 1; i < argc; i++) {
		CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
		CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
		CHECK_ARGUMENT_INT(i, "-n", numThreads , numThreadsSet)
	}

	if( !graphFileNameSet ) {
		printf("Graph filename not set\n");
		return 1;
	}

	if( !outputFileNameSet ) {
		outputFileName = new char_t[512];
		sprintf(outputFileName,"communities.dat");
	}

        CGraph graph;
	uint64_t totalTime = 0;
        uint64_t initTime, spendTime;

        
        //==================== LOAD THE GRAPH ==================================
	initTime = startClock();
	printf( "Graph: %s\n", graphFileName ); 
	printf( "OutputFile: %s\n", outputFileName );	
	graph.Load(graphFileName,numThreads);
	spendTime = stopClock(initTime);
	totalTime += spendTime;
	printf("Load time: %lu ms\n", spendTime);
        //======================================================================

        
	//================ REMOVE EDGES WITHOUT TRIANGLES ======================
	initTime = startClock();
	printf("Removing edges without triangles ...\n");
	graph.RemoveEdgesNoTriangles(numThreads);
        spendTime = stopClock(initTime);
	totalTime += spendTime;
	printf("Removing edges without triangles time: %lu ms\n", spendTime);
        //======================================================================

        
	//=================== INITIALIZE PARTITIONS ============================
	initTime = startClock();
	CommunityPartition partition;
	printf("Computing initial partition ...\n");
	if(InitializeSimplePartition(&graph,&partition)) {
		printf("Error computing initial partition\n");
		return 1;
	}
        spendTime = stopClock(initTime);
	totalTime += spendTime;
	printf("Initial partition time: %lu ms\n", spendTime);
        //======================================================================

        
        //================ TRANSFER NODES AMONG PARTITIONS =====================
        initTime = startClock();
        if(ImproveCommunities( &graph, &partition, numThreads) ) {
		printf("Error while improving communities\n");
		return 1;
	}
	spendTime = stopClock(initTime);
	totalTime += spendTime;
	printf("Improvement execution time: %lu ms\n", spendTime);
        //======================================================================

        
        //======================== PRINT RESULTS ===============================
        initTime = startClock();
	PrintPartition( &graph, &partition, outputFileName );
	spendTime = stopClock(initTime);
	totalTime += spendTime;
	printf("Print partition time: %lu ms\n", spendTime);
        //======================================================================


	printf("\n");
	printf("\n");
	printf("*******************************************************\n");
	printf("%-32s %-10d\n","Number of Communities:", partition.m_NumCommunities);
	printf("%-32s %-10f\n","WCC:", partition.m_WCC / (float32_t)graph.GetNumNodes() );
	printf("%-32s %-10lu ms\n", "Total execution time:", totalTime);
	printf("*******************************************************\n");

	FreeResources(&partition);

	if( !outputFileNameSet ) {
		delete [] outputFileName;
	}
	return 0;

}



