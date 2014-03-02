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
    printf("Usage: wcc <flags>\n");
    printf("Availaible flags:\n");
    printf("\t-f [network file name] : Specifies the network file.\n");
    printf("\t-p [partition file name] : Specifies the partition file name.\n");
}


int main(int argc, char ** argv) {

    bool graphFileNameSet = false;
    bool partitionFileNameSet = false;
    bool numThreadsSet = false;
    char_t * graphFileName = NULL;
    char_t * partitionFileName = NULL;
    uint32_t numThreads = omp_get_num_procs();

    for (uint32_t i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-p", partitionFileName, partitionFileNameSet)
    }

    if (!graphFileNameSet) {
        printf("Graph filename not set\n");
        PrintUsage();
        return 1;
    }
    
    if (!partitionFileNameSet) {
        printf("Partition filename not set\n");
        PrintUsage();
    return 1;
    }

    CGraph graph;

    //==================== LOAD THE GRAPH ==================================
    printf("Graph: %s\n", graphFileName);
    graph.Load(graphFileName, numThreads);
    graph.RemoveEdgesNoTriangles(numThreads);
    //======================================================================


    //=================== LOAD PARTITION ============================
    printf("PartitionFile: %s\n", partitionFileName);
    CommunityPartition partition;
    LoadPartition(&graph,&partition,partitionFileName);
    //======================================================================

    printf("*******************************************************\n");
    printf("%-32s %-10d\n", "Number of Communities:", partition.m_NumCommunities);
    printf("%-32s %-10f\n", "WCC:", partition.m_WCC / (float32_t) graph.GetNumNodes());
    printf("*******************************************************\n");

    FreeResources(&partition);
    return 0;
}



