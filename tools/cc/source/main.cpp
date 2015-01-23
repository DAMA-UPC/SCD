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
#include <wcc/wcc.h>
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
    bool alphaSet = false;
    char_t * graphFileName = NULL;
    char_t * partitionFileName = NULL;
    uint32_t numThreads = omp_get_num_procs();
    double alpha = 1.0;

    for (uint32_t i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
    }

    if (!graphFileNameSet) {
        printf("Graph filename not set\n");
        PrintUsage();
        return 1;
    }
    
    CGraph graph;

    //==================== LOAD THE GRAPH ==================================
    printf("Graph: %s\n", graphFileName);
    graph.Load(graphFileName, numThreads);

    double64_t cc = 0.0;
    for( int i = 0; i < graph.GetNumNodes(); ++i ) {
        uint32_t degree = graph.GetDegree(i);
        uint32_t numTriangles = 0;
        const uint32_t* adjacencies = graph.GetNeighbors(i);
        for( int j = 0; j < degree; ++j) {
            uint32_t neighbor = adjacencies[j];
            const uint32_t* adjacencies2 = graph.GetNeighbors(neighbor);
            uint32_t degree2 = graph.GetDegree(neighbor);
            uint32_t intersect = Intersect((uint32_t*)adjacencies,degree,(uint32_t*)adjacencies2,degree2);
            numTriangles+=intersect;
        }
        cc+= degree > 1 ? numTriangles / (double64_t)(degree*(degree-1)) : 0;
    }
    cc /= graph.GetNumNodes();
    printf("CC: %f\n", cc);
    //======================================================================
    return 0;
}



