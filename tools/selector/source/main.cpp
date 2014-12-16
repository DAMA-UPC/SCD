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
#include <wcc/wcc.h>
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
    printf("Usage: selector <flags>\n");
    printf("Availaible flags:\n");
    printf("\t-f [network file name] : Specifies the network file.\n");
    printf("\t-o [output file name] : Specifies the output file.\n");
    printf("\t-p [partition file name] : Specifies the partition file name.\n");
    printf("\t-a [alpha value] : Specifies the alpha value.\n");
    printf("\t-min_size [minimum size] : Specifies the minimum size of the communities to find.\n");
    printf("\t-max_size [max size] : Specifies the maximum size of the communities to find.\n");
    printf("\t-min_wcc [min wcc] : Specifies the minimum wcc of the communities to find.\n");
    printf("\t-max_wcc [max wcc] : Specifies the maximum wcc of the communities to find.\n");
}


int main(int argc, char ** argv) {

    bool graphFileNameSet = false;
    bool partitionFileNameSet = false;
    bool outputFileNameSet = false;
    bool numThreadsSet = false;
    bool alphaSet = false;
    bool minSizeSet = false;
    bool maxSizeSet = false;
    bool minWCCSet = false;
    bool maxWCCSet = false;
    uint32_t minSize= 1; 
    uint32_t maxSize= 10; 
    double  minWCC = 0.0; 
    double  maxWCC = 1.0; 
    double  alpha = 1.0;
    char_t * graphFileName = NULL;
    char_t * partitionFileName = NULL;
    char_t * outputFileName = NULL;
    uint32_t numThreads = omp_get_num_procs();

    for (uint32_t i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-p", partitionFileName, partitionFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
        CHECK_ARGUMENT_FLOAT(i, "-a", alpha, alphaSet)
        CHECK_ARGUMENT_INT(i, "-min_size", minSize, minSizeSet)
        CHECK_ARGUMENT_INT(i, "-max_size", maxSize, maxSizeSet)
        CHECK_ARGUMENT_FLOAT(i, "-min_wcc", minWCC, minWCCSet)
        CHECK_ARGUMENT_FLOAT(i, "-max_wcc", maxWCC, maxWCCSet)
    }

    if (!graphFileNameSet || !minSize || !maxSize || !minWCC || !maxWCC || !outputFileNameSet) {
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

    //=================== LOAD PARTITION ===================================
    printf("PartitionFile: %s\n", partitionFileName);
    CommunityPartition partition;
    LoadPartition(&graph,&partition,partitionFileName, alpha);
    //======================================================================

    std::cout << alpha << " " << minSize << " " << maxSize << " " << minWCC << " " << maxWCC << std::endl;

    std::ofstream outputFile;
    outputFile.open(outputFileName);

    for( int i = 0; i < partition.m_NumCommunities; ++i ) {
        std::set<uint32_t> community;
        uint32_t communitySize = partition.m_Communities[partition.m_CommunityIndices[i]];
        if( communitySize < minSize || communitySize > maxSize ) continue;
        uint32_t * community_ptr = &partition.m_Communities[partition.m_CommunityIndices[i]+1];
        for( int j = 0; j < communitySize; ++j ) {
            community.insert(community_ptr[j]);
        } 
        double score = ComputeWCC(&graph,alpha,community);
        if( score < minWCC || score > maxWCC ) continue;
        for( std::set<uint32_t>::iterator it = community.begin(); it != community.end(); ++it ) {
            outputFile << graph.ReMap(*it) << " ";
        }
        outputFile << std::endl;
        std::cout << score << std::endl;
    }

    outputFile.close();
    FreeResources(&partition);
    return 0;
}
