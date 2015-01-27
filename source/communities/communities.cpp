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

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <common/time.h>
#include <communities/communities.h>
#include <fstream>
#include <graph/graph.h>
#include <iostream>
#include <map>
#include <omp.h>
#include <sstream>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <wcc/wcc.h>

namespace scd {

#define SCD_INVALID_COMMUNITY 0xffffffff

    uint32_t num_threads = 1;

    /** @brief Types of movements.*/
    enum MovementType {
        E_REMOVE,
        E_REMOVE_AND_INSERT,
        E_NO_MOVEMENT
    };

    /** @brief This struct represents a movement.*/
    struct Movement {
        MovementType m_MovementType;
        uint32_t m_NodeId;
        uint32_t m_Community;
        double64_t m_Improvement;
    };

    /** @brief Compares two node clustering.
     *  @param e1 Void pointer to the first node clustering.    
     *  @param e2 Void pointer to the second node clustering.   
     *  @return -1 if e1 goes before e2. 1 if e1 goes after e2. 0 if e1 and e2 are equal.*/
    static int Compare_NodeClusterings(const void* e1, const void* e2) {
        NodeClustering* nC1 = (NodeClustering*) e1;
        NodeClustering* nC2 = (NodeClustering*) e2;
        if (nC1->m_CC > nC2->m_CC) return -1;
        if (nC1->m_CC < nC2->m_CC) return 1;
        if (nC1->m_Degree > nC2->m_Degree) return -1;
        if (nC1->m_Degree < nC2->m_Degree) return 1;
        return 0;
    }

    /** @brief Compares two unsigned integers.
     *  @param e1 Void pointer to the first unsigned integer.   
     *  @param e2 Void pointer to the second unsigned integer.  
     *  @return -1 if e1 goes before e2. 1 if e1 goes after e2. 0 if e1 and e2 are equal.*/
    static int Compare_Ids(const void* e1, const void* e2) {
        uint32_t id1 = *(uint32_t*) e1;
        uint32_t id2 = *(uint32_t*) e2;
        if (id1 < id2) return -1;
        if (id2 < id1) return 1;
        return 0;
    }

    /** @brief this function is used to compress the laberl space of a partition in order to be comprissed between
     *           0 and the actual number of communities - 1.
     *  @param[in] graph The graph where the partition belongs.
     *  @param[in] communities The array of current labels.
     *  @param[out] destCommunities The array where the new labeling will be stored.*/
    static uint32_t CompressCommunityLabels(const CGraph* graph, const uint32_t * communities, uint32_t * destCommunities) {
        std::map<uint32_t, uint32_t> * map = new std::map<uint32_t, uint32_t>();
        uint32_t label = 0;
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            if (!map->count(communities[i])) {
                map->insert(std::pair<uint32_t, uint32_t>(communities[i], label));
                label++;
            }
        }

        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            destCommunities[i] = (*(map->find(communities[i]))).second;
        }
        delete map;
        return label;
    }

    /** @brief Initializes a partition structure from a labels to communities array.
     *  @param[in] graph The graph where the partition belongs.
     *  @param[out] partition The partition structure where the partition will be stored.
     *  @param[in] communities The array of nodes to community labels from which the partition is initialized.
     *  @param[in] alfa The alfa parameter controlling the cohesivness of the communities.*/
    static uint32_t InitializeFromLabelsArray(const CGraph* graph, CommunityPartition* partition, const uint32_t* communities, const double64_t alfa) {

        //Initializing default values
        partition->m_NodeLabels = NULL;
        partition->m_CommunityIndices = NULL;
        partition->m_Communities = NULL;
        partition->m_InternalEdges = NULL;
        partition->m_ExternalEdges = NULL;
        partition->m_NodeWCC = NULL;
        partition->m_NumCommunities = 0;
        partition->m_WCC = 0;

        partition->m_NumNodes = graph->GetNumNodes();
        partition->m_NodeLabels = new uint32_t[graph->GetNumNodes()];
        if (!partition->m_NodeLabels) {
            printf("Error allocating node labels %u.\n", graph->GetNumNodes());
            return 1;
        }

        uint32_t maxNumCommunities = 0;
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            if (communities[i] > maxNumCommunities) {
                maxNumCommunities = communities[i];
            }
        }
        maxNumCommunities++;

        partition->m_NumCommunities = CompressCommunityLabels(graph, communities, partition->m_NodeLabels);


        //Allocating space to store the communities
        partition->m_CommunityIndices = new uint32_t[partition->m_NumCommunities];
        if (!partition->m_CommunityIndices) {
            printf("Error allocating labels indices.\n");
            return 1;
        }
        partition->m_Communities = new uint32_t[partition->m_NumCommunities + graph->GetNumNodes()];
        if (!partition->m_Communities) {
            printf("Error allocating inverted index.\n");
            return 1;
        }

        partition->m_NodeWCC = new double64_t[graph->GetNumNodes()];
        if (!partition->m_NodeWCC) {
            printf("Error allocating node labels %u.\n", graph->GetNumNodes());
            return 1;
        }

        //Creating the counters the creation of the inverted index
        uint32_t* counters = new uint32_t[partition->m_NumCommunities];
        if (!counters) {
            printf("Error allocating counters: %u\n", partition->m_NumCommunities);
            return 1;
        }

#pragma omp parallel for schedule(SCD_SCHEDULING, SCD_THREAD_BLOCK_SIZE)
        for (uint32_t i = 0; i < partition->m_NumCommunities; i++) {
            counters[i] = 0;
        }
        //Computing community sizes;
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            counters[partition->m_NodeLabels[i]]++;
        }
        //Initializing labels indices.
        uint32_t currentIndex = 0;
        for (uint32_t i = 0; i < partition->m_NumCommunities; i++) {
            if (counters[i] > 0) {
                partition->m_CommunityIndices[i] = currentIndex;
                partition->m_Communities[currentIndex] = counters[i];
                currentIndex += counters[i] + 1;
            } else {
                partition->m_CommunityIndices[i] = SCD_INVALID_COMMUNITY;
            }
            counters[i] = 0;
        }
        //Initializing the inverted index.
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            uint32_t lIndex = partition->m_CommunityIndices[partition->m_NodeLabels[i]];
            assert(lIndex != SCD_INVALID_COMMUNITY);
            assert(counters[partition->m_NodeLabels[i]] < partition->m_Communities[lIndex]);
            partition->m_Communities[lIndex + counters[partition->m_NodeLabels[i]] + 1] = i;
            counters[partition->m_NodeLabels[i]]++;
        }

        for (uint32_t i = 0; i < partition->m_NumCommunities; i++) {
            if (partition->m_CommunityIndices[i] != SCD_INVALID_COMMUNITY) {
                uint32_t lIndex = partition->m_CommunityIndices[i];
                qsort(&(partition->m_Communities[lIndex + 1]), partition->m_Communities[lIndex], sizeof (uint32_t), Compare_Ids);
            }
        }
        delete[] counters;

        partition->m_InternalEdges = new uint32_t[partition->m_NumCommunities];
        if (!partition->m_InternalEdges) {
            printf("Error while allocating internal edges.\n");
            return 1;
        }

        partition->m_ExternalEdges = new uint32_t[partition->m_NumCommunities];
        if (!partition->m_ExternalEdges) {
            printf("Error while allocating external edges.\n");
            return 1;
        }

#pragma omp parallel for schedule(SCD_SCHEDULING, SCD_THREAD_BLOCK_SIZE)
        for (uint32_t i = 0; i < partition->m_NumCommunities; i++) {
            partition->m_InternalEdges[i] = 0;
            partition->m_ExternalEdges[i] = 0;
        }

        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            const uint32_t* adjacencies = graph->GetNeighbors(i);
            uint32_t degree = graph->GetDegree(i);
            for (uint32_t j = 0; j < degree; j++) {
                if (i < adjacencies[j]) {
                    if (partition->m_NodeLabels[i] == partition->m_NodeLabels[adjacencies[j]]) {
                        partition->m_InternalEdges[partition->m_NodeLabels[i]]++;
                    } else {
                        partition->m_ExternalEdges[partition->m_NodeLabels[i]]++;
                        partition->m_ExternalEdges[partition->m_NodeLabels[adjacencies[j]]]++;
                    }
                }
            }
        }

        partition->m_WCC = ComputeWCC(graph, alfa, partition->m_NodeLabels, partition->m_NumCommunities, partition->m_CommunityIndices, partition->m_Communities, partition->m_NodeWCC);
        return 0;
    }

    /** @brief Computes the increment on WCC for inserting a node into a community.
      @param[in] r The size of the community.
      @param[in] d_in The number of edges between the inserted vertex and the community.
      @param[in] d_out The number of edges between the inserted vertex and the rest of the graph.
      @param[in] c_out The number of edges leaving the community (note that this MUST include d_in).
      @param[in] p_in The probability that an edge inside of the community exists.
      @param[in] p_ext The probability that two edges leaving the community close a triangle.
      @param[in] alfa The alfa parameter controlling the cohesivness of the communities.*/
    static double64_t CheckForIncrement(int32_t r, int32_t d_in, int32_t d_out, uint32_t c_out, double64_t p_in, double64_t p_ext, const double64_t alfa) {
       /* double64_t t;
        if (r > 0) {
            t = (c_out - d_in) / (double64_t) r;
        } else {
            t = 0.0;
        }
        double64_t A = 0.0;
        double64_t denom = 0.0;
        denom = (d_in * (d_in - 1) * p_in + d_out * (d_out + d_in - 1) * p_ext);
        if (denom != 0.0 && ((r + d_out) > 0)) {
            A = ((d_in * (d_in - 1) * p_in) / denom) * (d_in + d_out) / (double64_t) (r + d_out);
        }
        double64_t BMinus = 0.0;
        denom = (r - 1)*(r - 2) * p_in * p_in * p_in + (d_in - 1) * p_in + t * (r - 1) * p_in * p_ext + t * (t - 1) * p_ext + (d_out) * p_ext;
        if (denom != 0.0 && ((r + t) > 0)) {
            BMinus = (((d_in - 1) * p_in) / denom) * ((r - 1) * p_in + 1 + t) / (r + t);
        }
        double64_t CMinus = 0.0;
        denom = (r - 1)*(r - 2) * p_in * p_in * p_in + t * (t - 1) * p_ext + t * (r - 1)*(p_in) * p_ext;
        if (denom != 0.0 && ((r + t) > 0) && ((r - 1 + t) > 0)) {
            CMinus = -(((r - 1)*(r - 2) * p_in * p_in * p_in) / denom) * ((r - 1) * p_in + t) / ((r + t)*(r - 1 + t));
        }
        return (A + d_in * BMinus + (r - d_in) * CMinus);
        */

        double64_t t;
        if (r > 0) {
            t = (c_out - d_in) / (double64_t) r;
        } else {
            t = 0.0;
        }
        // Node v 
        double64_t A = 0.0;
        double64_t denom = 0.0;
        denom = (d_in * (d_in - 1) * p_in + 
                d_out * (d_out - 1) * p_ext) +
                d_out * d_in * p_ext;
        denom *= d_in + d_out + alfa*(r-1-d_in);
        if (denom != 0.0) {
            A = ((d_in * (d_in - 1) * p_in) * (d_in + d_out)) / (denom);
        }

        // Nodes connected with v 
        double64_t BMinus = 0.0;
        denom = (r - 1)*(r - 2) * p_in * p_in * p_in + 
                2*(d_in - 1) * p_in + 
                t * (r - 1) * p_in * p_ext + 
                t * (t - 1) * p_ext + 
                (d_out) * p_ext;
        denom *= (r-1)*p_in + 1 + t + alfa*(r - (r-1)*p_in - 1);
        if (denom != 0.0) {
            BMinus = (2*(d_in - 1) * p_in) * ((r - 1) * p_in + 1 + t) / denom;
        }
        // Nodes not connected with v 
        double64_t CMinus = 0.0;
        denom = (r - 1)*(r - 2) * p_in * p_in * p_in +
                t * (t - 1) * p_ext + 
                t * (r - 1)*(p_in) * p_ext;
        denom *= (r-1)*p_in + t + alfa*(r - (r-1)*p_in);
        denom *= (r-1)*p_in + t + alfa*(r - (r-1)*p_in - 1);
        if (denom != 0.0 && ((r + t) > 0) && ((r - 1 + t) > 0)) {
            CMinus = -((r - 1)*(r - 2) * p_in * p_in * p_in) * ((r - 1) * p_in + t)*alfa / denom;
        }
        // Total 
        return (A + d_in * BMinus + (r - d_in) * CMinus);
    }

    /** @brief Checks the best movement of a vertex.
      @param[in] graph The graph.
      @param[in] node The node to check the movement.
      @param[in] partition The current partition into communities.
      @return The movement to perform.*/
    static Movement CheckForBestMovement(const CGraph* graph, uint32_t node, const CommunityPartition* partition, const double64_t alfa) {

        Movement movement;
        movement.m_MovementType = E_NO_MOVEMENT;
        movement.m_NodeId       = node;

        std::map<uint32_t, uint32_t> neighborsCommunity;
        neighborsCommunity.insert(std::pair<uint32_t, uint32_t>(partition->m_NodeLabels[node], 0));
        const uint32_t * adjacencies = graph->GetNeighbors(node);
        uint32_t degree = graph->GetDegree(node);
        for (uint32_t i = 0; i < degree; i++) {
            uint32_t neighbor = adjacencies[i];
            if (partition->m_Communities[partition->m_CommunityIndices[partition->m_NodeLabels[neighbor]]] > 1) {
                std::map<uint32_t, uint32_t>::iterator it = neighborsCommunity.find(partition->m_NodeLabels[neighbor]);
                if (it != neighborsCommunity.end()) {
                    (*it).second++;
                } else {
                    neighborsCommunity.insert(std::pair<uint32_t, uint32_t>(partition->m_NodeLabels[neighbor], 1));
                }
            }
        }

        bool removeCommunity = false;
        uint32_t bestRemoveInternalEdges = 0;
        uint32_t auxInternalEdges = (*neighborsCommunity.find(partition->m_NodeLabels[node])).second;
        uint32_t community = partition->m_NodeLabels[node];
        uint32_t communityIndex = partition->m_CommunityIndices[community];
        uint32_t communitySize = partition->m_Communities[communityIndex];
        double64_t p_in;
        if ((communitySize - 2) != 0 && (communitySize - 1) != 0) {
            p_in = (2 * partition->m_InternalEdges[community] - auxInternalEdges * 2) / ((double64_t) (communitySize - 1) * (communitySize - 2));
        } else {
            p_in = 0.0f;
        }
        double64_t p_ext = graph->GetCC();
        double64_t bestRemoveImprovement;
        bestRemoveImprovement = -CheckForIncrement(communitySize - 1, auxInternalEdges,
                graph->GetDegree(node) - auxInternalEdges,
                partition->m_ExternalEdges[community] + auxInternalEdges - (graph->GetDegree(node) - auxInternalEdges ),
                p_in, p_ext, alfa);
        bestRemoveInternalEdges = auxInternalEdges;
        if (bestRemoveImprovement > 0.0f) {
            removeCommunity = true;
        }


        uint32_t   bestInsertCommunity;
        double64_t bestInsertImprovement      = -10000000000000.0;
        uint32_t   bestInsertInternalEdges    = 0;
        bool       insertCommunity            = false;

        for (std::map<uint32_t, uint32_t>::iterator it = neighborsCommunity.begin(); it != neighborsCommunity.end(); ++it) {
            uint32_t community = it->first;
            if (community != partition->m_NodeLabels[node]) {
                uint32_t auxInternalEdges = it->second;
                uint32_t communityIndex   = partition->m_CommunityIndices[community];
                uint32_t communitySize    = partition->m_Communities[communityIndex];
                double64_t p_in;
                if ((communitySize - 1) > 0 && (communitySize > 0)) {
                    p_in = (2 * partition->m_InternalEdges[community]) / ((double64_t) (communitySize) * (communitySize - 1));
                } else {
                    p_in = 0.0;
                }

                double64_t p_ext          = graph->GetCC();
                double64_t auxImprovement = CheckForIncrement(communitySize, auxInternalEdges, graph->GetDegree(node) - auxInternalEdges,
                        partition->m_ExternalEdges[community], p_in, p_ext, alfa);
                if (auxImprovement + bestRemoveImprovement > bestInsertImprovement) {
                    insertCommunity = true;
                    bestInsertImprovement = auxImprovement + bestRemoveImprovement;
                    bestInsertCommunity = community;
                    bestInsertInternalEdges = auxInternalEdges;
                }
            }
        }

        if (bestInsertImprovement > 0.0f && ((bestInsertImprovement) > bestRemoveImprovement)) {
            movement.m_MovementType = E_REMOVE_AND_INSERT;
            movement.m_Community = bestInsertCommunity;
            movement.m_Improvement = bestInsertImprovement;
        } else if (bestRemoveImprovement > 0.0f) {
            movement.m_MovementType = E_REMOVE;
            movement.m_Improvement = bestRemoveImprovement;
        }
        return movement;
    }


    /** @brief Performs an improvement step, that is, checks for movements for all the nodes and
      and computes the new partitions.
      @param[in] graph The graph.
      @param[out] partition The current partition. It will be modified with the new partition.
      @param[in] alfa The alfa parameter controlling the cohesivness of the communities.*/
    static uint32_t PerformImprovementStep(const CGraph* graph, CommunityPartition* partition, const double64_t alfa) {
        std::vector<Movement>* movements = new std::vector<Movement>[num_threads];
        uint32_t N = graph->GetNumNodes();

#pragma omp parallel for schedule(SCD_SCHEDULING,SCD_THREAD_BLOCK_SIZE) 
        for (uint32_t i = 0; i < N; i++) {
            int thread = omp_get_thread_num();
            if (i % 100000 == 0) {
                printf("Thread %d: Checked movements of %d nodes.\n", thread, i);
            }
            Movement movement;
            movement = CheckForBestMovement(graph, i, partition, alfa);
            if (movement.m_MovementType != E_NO_MOVEMENT) {
                movements[thread].push_back(movement);
            }
        }
        printf("All movements checked\n");

        uint32_t* tempNodeLabels = new uint32_t[partition->m_NumNodes];
        memcpy(&tempNodeLabels[0], &partition->m_NodeLabels[0], sizeof (uint32_t) * partition->m_NumNodes);
        uint32_t totalMovements = 0;

        //uint32_t nextLabel = partition->m_NumCommunities;
        uint32_t removeMovements = 0;
        uint32_t removeAndInsertMovements = 0;
        uint32_t insertMovements = 0;


#pragma omp parallel for schedule(static,1)   
        for (uint32_t thread = 0; thread < num_threads; thread++) {
            uint32_t numMovements = movements[thread].size();
            totalMovements += numMovements;
            uint32_t nextLabelThread = partition->m_NumCommunities + numMovements * thread;            

            for (uint32_t i = 0; i < numMovements; i++) {
                Movement movement = (movements[thread])[i];
                switch (movement.m_MovementType) {
                    case E_REMOVE:
                        tempNodeLabels[movement.m_NodeId] = nextLabelThread;
                        removeMovements++;
                        nextLabelThread++;
                        break;
                    case E_REMOVE_AND_INSERT:
                        tempNodeLabels[movement.m_NodeId] = movement.m_Community;
                        if (partition->m_Communities[partition->m_CommunityIndices[partition->m_NodeLabels[movement.m_NodeId]]] == 1) {
                            insertMovements++;
                        } else {
                            removeAndInsertMovements++;
                        }
                        break;
                }
            }
        }
        delete [] movements;
        printf(" Number of removes performed: %d\n", removeMovements);
        printf(" Number of remove and insert performed: %d\n", removeAndInsertMovements);
        printf(" Number of insert performed: %d\n", insertMovements);
        FreeResources(partition);

        if (InitializeFromLabelsArray(graph, partition, tempNodeLabels, alfa)) {
            printf("Error initializing from label array.\n");
            return 1;
        }
        delete [] tempNodeLabels;

        return 0;
    }


    /*********************** EXPERIMENTAL **************************************/

    struct CommunityInteraction {
        uint32_t m_CommunityId1;
        uint32_t m_CommunityId2;
        uint32_t degree;
        double64_t m_Improvement;
    };

    static bool CompareByImprovement( const CommunityInteraction& a, const CommunityInteraction& b ) {
        if( a.m_Improvement > b.m_Improvement ) return true;
        return false;
    }

    static bool CompareById( const CommunityInteraction& a, const CommunityInteraction& b ) {
        if( a.m_CommunityId1 < b.m_CommunityId1 ) return true;
        if( a.m_CommunityId1 > b.m_CommunityId1 ) return false;
        if( a.m_CommunityId2 < b.m_CommunityId2 ) return true;
        if( a.m_CommunityId2 > b.m_CommunityId2 ) return false;
        return false;
    }

    static void PrintCommunity(const CommunityPartition* partition, uint32_t communityId) {
        const uint32_t* nodes = &partition->m_Communities[partition->m_CommunityIndices[communityId]+1];
        uint32_t size = partition->m_Communities[partition->m_CommunityIndices[communityId]];
        for( uint32_t i = 0; i < size; ++i ) {
            printf("%d ", nodes[i]);
        }
        printf("\n");
    }

    double64_t TestMerge(const CGraph* graph, const CommunityPartition* partition, const double64_t alfa, const CommunityInteraction& interaction) {
//        printf("Testing merge %d - %d\n", interaction.m_CommunityId1, interaction.m_CommunityId2);
//        PrintCommunity(partition, interaction.m_CommunityId1);
//        PrintCommunity(partition, interaction.m_CommunityId2);

        std::set<uint32_t> community;
        const uint32_t* nodes = &partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId1]+1];
        uint32_t size = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId1]];
        double64_t before = 0.0;
        for( uint32_t i = 0; i < size; ++i ) {
            community.insert(nodes[i]);
            before += partition->m_NodeWCC[nodes[i]];
        }

        nodes = &partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId2]+1];
        size = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId2]];
        for( uint32_t i = 0; i < size; ++i ) {
            community.insert(nodes[i]);
            before += partition->m_NodeWCC[nodes[i]];
        }
 /*       for( std::set<uint32_t>::iterator it = community.begin(); it != community.end(); ++it){
            printf("%d ",*it);
        }
        printf("\n");*/
        double64_t after = ComputeWCC(graph,alfa,community);
//        printf("after: %f, before: %f\n", after, before);
        return (after - before);
    }


    void ComputeDegree( const CGraph* graph, const CommunityPartition* partition, const double64_t alfa, CommunityInteraction& interaction) {
        interaction.degree = 0;
        uint32_t communitySize1 = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId1]];
        uint32_t communitySize2 = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId2]];
        uint32_t* community1 = &partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId1]+1];
        uint32_t* community2 = &partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId2]+1];
        for(int i = 0; i < communitySize1; ++i) {
            uint32_t node = community1[i]; 
            const uint32_t* adjacencies = graph->GetNeighbors(node); 
            uint32_t degree = graph->GetDegree(node);
            for( int j = 0; j < degree;++j){
                uint32_t neighbor = adjacencies[j];
                if( partition->m_NodeLabels[neighbor] == interaction.m_CommunityId2 ) {
                    interaction.degree++;
                }
            }
        }
    }

    double ShouldMerge( const CGraph* graph, const CommunityPartition* partition, const double64_t alfa, CommunityInteraction& interaction) {
        ComputeDegree(graph,partition,alfa,interaction);
        uint32_t communitySize1 = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId1]];
        uint32_t communitySize2 = partition->m_Communities[partition->m_CommunityIndices[interaction.m_CommunityId2]];
        double64_t p_in1 = (2 * partition->m_InternalEdges[interaction.m_CommunityId1]) / ((double64_t) (communitySize1) * (communitySize1 - 1));
        double64_t p_in2 = (2 * partition->m_InternalEdges[interaction.m_CommunityId2]) / ((double64_t) (communitySize2) * (communitySize2 - 1));
        double64_t p_in = (p_in1*communitySize1 + p_in2*communitySize2) / (communitySize1+communitySize2);
        //double64_t p_in = std::max(p_in1,p_in2);
        double64_t p_out = interaction.degree / (double64_t)(communitySize1*communitySize2);
        double64_t threshold = - (-alfa*p_out + sqrt(alfa*alfa*p_out*p_out-2*alfa*alfa*p_out-4*alfa*p_out*p_out+2*alfa*alfa+2*alfa*p_out+p_out*p_out)+p_out)*p_out / alfa*(p_out-1.0); 
        //std::cout << p_in << " " << p_out << " " << threshold << std::endl;
        return threshold - p_in;
    }

    typedef std::set<CommunityInteraction, bool (*)(const CommunityInteraction&, const CommunityInteraction&)> InteractionsSet;
    static uint32_t MergeCommunities(const CGraph* graph, CommunityPartition* partition, const double64_t alfa) {
        // Look for community interactions.
        InteractionsSet candidateMerges(CompareById);
        uint32_t N = graph->GetNumNodes();
        for( uint32_t i = 0; i < N; ++i ) {
            uint32_t communityLabel1 = partition->m_NodeLabels[i];
            uint32_t degree = graph->GetDegree(i);
            const uint32_t* adjacencies = graph->GetNeighbors(i);
            for( uint32_t j = 0; j < degree; ++j ) {
                uint32_t communityLabel2 = partition->m_NodeLabels[adjacencies[j]];
                if( communityLabel1 != communityLabel2 ) {
                    CommunityInteraction cI;
                    cI.degree = 0;
                    if( communityLabel1 < communityLabel2 ) {
                        cI.m_CommunityId1 = communityLabel1;
                        cI.m_CommunityId2 = communityLabel2;
                    } else {
                        cI.m_CommunityId2 = communityLabel1;
                        cI.m_CommunityId1 = communityLabel2;
                    }
                    candidateMerges.insert(cI);
                }
            }
        }
        std::vector<CommunityInteraction> filteredInteractions;
        uint32_t earlyFilter = 0;
        // Test each community interaction and rank it.
        for( InteractionsSet::iterator it = candidateMerges.begin(); it != candidateMerges.end(); ++it ) {
            CommunityInteraction cI = *it;
            uint32_t size1 = partition->m_Communities[partition->m_CommunityIndices[cI.m_CommunityId1]];
            uint32_t size2 = partition->m_Communities[partition->m_CommunityIndices[cI.m_CommunityId2]];
            /*if( size1 == 1 || size2 == 1 )*/ {
                /*if(ShouldMerge(graph, partition, alfa, cI) > 0.0 )*/ {
                    //ShouldMerge(graph, partition, alfa, cI);
                    earlyFilter++;
                    double64_t improvement = TestMerge(graph, partition, alfa, cI);
                    if( improvement > 0.0 ) {
//                        std::cout << size1 << " " << size2 << " " << cI.degree << std::endl;
                        //std::cout << cI.m_CommunityId1 << " " << cI.m_CommunityId2 << " " << improvement << std::endl;
                        cI.m_Improvement = improvement;
                        filteredInteractions.push_back(cI);
                    }
                }
            }
        } 
        std::cout << earlyFilter << " " << filteredInteractions.size() << " " << candidateMerges.size() << std::endl;
        // Sort community interactions by improvement.
        uint32_t* tempNodeLabels = new uint32_t[partition->m_NumNodes];
        memcpy(tempNodeLabels, partition->m_NodeLabels, sizeof(uint32_t)*partition->m_NumNodes);
        std::sort(filteredInteractions.begin(), filteredInteractions.end(), CompareByImprovement);
        std::set<uint32_t> touched;
        uint32_t numInteractions = filteredInteractions.size();
        for( uint32_t i = 0; i < numInteractions; ++i ) {
            if( (touched.find(filteredInteractions[i].m_CommunityId1) == touched.end()) && 
                    (touched.find(filteredInteractions[i].m_CommunityId2) == touched.end()) ) {
                uint32_t communitySize = partition->m_Communities[partition->m_CommunityIndices[filteredInteractions[i].m_CommunityId1]];
                const uint32_t* community = &partition->m_Communities[partition->m_CommunityIndices[filteredInteractions[i].m_CommunityId1]+1];
                for( uint32_t j = 0; j < communitySize; ++j ) {
                    tempNodeLabels[community[j]] = filteredInteractions[i].m_CommunityId2;
                }
                touched.insert(filteredInteractions[i].m_CommunityId1);
                touched.insert(filteredInteractions[i].m_CommunityId2);
            }
        }

        // Perform interactions constrained by independence and create a new labels array to create a partition from.
        FreeResources(partition);
        if (InitializeFromLabelsArray(graph, partition, tempNodeLabels, alfa)) {
            printf("Error initializing from label array.\n");
            return 1;
        }
        delete [] tempNodeLabels;
        return 0;
    }
    /***************************************************************************/

    /*****************************EXPERIMENTAL 2 ******************************/

    bool CompareMovements(const Movement& a, const Movement& b) {
        if(a.m_Community < b.m_Community) return true;
        if(b.m_Community < a.m_Community) return false;
        if(a.m_Improvement > b.m_Improvement) return true;
        if(b.m_Improvement > a.m_Improvement) return false;
        return false;
    }

    static uint32_t PerformImprovementStep2(const CGraph* graph, CommunityPartition* partition, const double64_t alfa) {
        std::vector<Movement>* movements = new std::vector<Movement>[num_threads];
        uint32_t N = graph->GetNumNodes();

#pragma omp parallel for schedule(SCD_SCHEDULING,SCD_THREAD_BLOCK_SIZE) 
        for (uint32_t i = 0; i < N; i++) {
            int thread = omp_get_thread_num();
            if (i % 100000 == 0) {
                printf("Thread %d: Checked movements of %d nodes.\n", thread, i);
            }
            Movement movement;
            movement = CheckForBestMovement(graph, i, partition, alfa);
            if (movement.m_MovementType != E_NO_MOVEMENT) {
                movements[thread].push_back(movement);
            }
        }
        printf("All movements checked\n");

        for( uint32_t i = 0; i < N; i++) {
            std::sort((movements[i]).begin(), (movements[i]).end(),CompareMovements);
        }

        uint32_t* tempNodeLabels = new uint32_t[partition->m_NumNodes];
        memcpy(&tempNodeLabels[0], &partition->m_NodeLabels[0], sizeof (uint32_t) * partition->m_NumNodes);
        uint32_t totalMovements = 0;

        //uint32_t nextLabel = partition->m_NumCommunities;
        uint32_t removeMovements = 0;
        uint32_t removeAndInsertMovements = 0;
        uint32_t insertMovements = 0;

        
#pragma omp parallel for schedule(static,1)   
        for (uint32_t thread = 0; thread < num_threads; thread++) {
            uint32_t numMovements = movements[thread].size();
            totalMovements += numMovements;
            uint32_t nextLabelThread = partition->m_NumCommunities + numMovements * thread;            

            uint32_t previousCommunity = 100000000;
            for (uint32_t i = 0; i < numMovements; i++) {
                Movement movement = (movements[thread])[i];
                if(movement.m_Community != previousCommunity) {
                    previousCommunity = movement.m_Community;
                    switch (movement.m_MovementType) {
                        case E_REMOVE:
                            tempNodeLabels[movement.m_NodeId] = nextLabelThread;
                            removeMovements++;
                            nextLabelThread++;
                            break;
                        case E_REMOVE_AND_INSERT:
                            tempNodeLabels[movement.m_NodeId] = movement.m_Community;
                            if (partition->m_Communities[partition->m_CommunityIndices[partition->m_NodeLabels[movement.m_NodeId]]] == 1) {
                                insertMovements++;
                            } else {
                                removeAndInsertMovements++;
                            }
                            break;
                    }
                }
            }
        }
        delete [] movements;
        printf(" Number of removes performed: %d\n", removeMovements);
        printf(" Number of remove and insert performed: %d\n", removeAndInsertMovements);
        printf(" Number of insert performed: %d\n", insertMovements);
        FreeResources(partition);

        if (InitializeFromLabelsArray(graph, partition, tempNodeLabels, alfa)) {
            printf("Error initializing from label array.\n");
            return 1;
        }
        delete [] tempNodeLabels;

        return 0;
    }

/*    Movement PerformCommunityCentricRefinement( const CGraph* graph, const CommunityPartition partition*, const double64_t alfa, uint32_t* communityId, uint32_t size ) {
        Movement bestMovement;
        bestMovement.m_Type = E_NO_MOVEMENT;
        for( int i = 0; i < size; ++i ) {
           uint32_t node = community[i];
           const uint32_t* adjacencies  = graph->GetNeighbors(node);
           const uint32_t degree = graph->GetDegree(node);
           for(uint32_t j = 0; j < degree; ++j ) {
               uint32_t neighbor = adjacencies[j];
               if(partition->m_Communities[partition->m_CommunityIndices[partition->m_NodeLabels[neighbor]]] == 1 ) {

               }
           }
        }
    }
    */



    /** @brief Measures the memory consumption of a partition.
      @param partition The partition to measure.
      @return The size in bytes of the structure.*/
    static uint64_t MeasureMemoryConsumption(const CommunityPartition* partition) {
        uint64_t memoryConsumption = 0;
        memoryConsumption += sizeof (uint32_t)  * partition->m_NumNodes; //Labels array consumption.
        memoryConsumption += sizeof (uint32_t)  * partition->m_NumCommunities; //Community indices consumption.
        memoryConsumption += sizeof (uint32_t)  *(partition->m_NumCommunities + partition->m_NumNodes); //Communities consumption.
        memoryConsumption += sizeof (uint32_t)  * partition->m_NumCommunities; //Internal edges consumption.
        memoryConsumption += sizeof (uint32_t)  * partition->m_NumCommunities; //External edges consumption.
        memoryConsumption += sizeof (double64_t)* partition->m_NumNodes; //WCCs consumption.
        memoryConsumption += sizeof (uint32_t); //NumNodes consumption.
        memoryConsumption += sizeof (uint32_t); //NumCommunities consumption.
        memoryConsumption += sizeof (double64_t); //WCC consumption.
        return memoryConsumption;
    }

    uint32_t    LoadPartition( const CGraph* graph, CommunityPartition* partition, const char_t* partitionFileName, const double64_t alfa ) {

        std::map<uint32_t, uint32_t> oldToNew;
        const uint32_t* newToOld = graph->GetMap();
        for( uint32_t i = 0; i < graph->GetNumNodes(); ++i ) {
            oldToNew.insert(std::pair<uint32_t, uint32_t>(newToOld[i],i));
        }

        uint32_t* communities = new uint32_t[graph->GetNumNodes()];
        memset(communities,0xff,sizeof(uint32_t)*graph->GetNumNodes());
        if (!communities) {
            printf("Unable to allocate partition\n");
            return 1;
        }

        std::ifstream partitionFile;
        partitionFile.open(partitionFileName);
        if(!partitionFile.is_open()) {
            printf("Unable to load partition file.\n");
            return 1;
        }
        std::string line;
        uint32_t nextLabel = 0;
        while(std::getline(partitionFile,line)) {
            std::istringstream stream(line);
            uint32_t node;
            while( stream >> node ) {
                communities[oldToNew[node]] = nextLabel;
            }
            ++nextLabel;
        }
        for( uint32_t i=0; i<graph->GetNumNodes(); ++i) {
            if( communities[i] == 0xffffffff ) {
                communities[i] = nextLabel++;
            }
        }
        partitionFile.close();
        InitializeFromLabelsArray(graph,partition,communities, alfa);
        delete [] communities;
        return 0;
    }

    uint32_t InitializeSimplePartition(const CGraph* graph, CommunityPartition* partition, const double64_t alfa) {
        //Computing the clustering coefficient of each node of the graph.
        NodeClustering* nC = new NodeClustering[graph->GetNumNodes()];
        if (!nC) {
            printf("Error allocating node clustering array.");
            return 1;
        }
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            nC[i].m_NodeId = i;
            nC[i].m_Degree = graph->GetDegree(i);
            nC[i].m_CC = graph->GetTotalTriangles(i) / (double64_t) (nC[i].m_Degree * (nC[i].m_Degree - 1));
        }
        qsort(nC, graph->GetNumNodes(), sizeof (NodeClustering), Compare_NodeClusterings);
        //Creating a vector to track which nodes have already been visited.
        bool * visited = new bool[graph->GetNumNodes()];
        if (!visited) {
            printf("Error allocating visited array.");
            return 1;
        }

        memset(visited, false, graph->GetNumNodes() );

        uint32_t* communities = new uint32_t[graph->GetNumNodes()];
        if (!communities) {
            printf("Error allocating communities array.\n");
            return 1;
        }

        uint32_t nextLabel = 0;
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            NodeClustering* nodeClustering = &nC[i];
            //          printf("%f\n", nodeClustering->m_CC);
            if (!visited[nodeClustering->m_NodeId]) {
                visited[nodeClustering->m_NodeId] = true;
                communities[nodeClustering->m_NodeId] = nextLabel;
                const uint32_t* adjacencies1 = graph->GetNeighbors(nodeClustering->m_NodeId);
                uint32_t degree = graph->GetDegree(nodeClustering->m_NodeId);

                for (uint32_t j = 0; j < degree; j++) {
                    if (!visited[adjacencies1[j]]) {                        
                        visited[adjacencies1[j]] = true;
                        communities[adjacencies1[j]] = nextLabel;                       
                    }
                }
                nextLabel++;
            }
        }
        delete [] visited;
        delete [] nC;

        InitializeFromLabelsArray(graph, partition, communities, alfa);
        delete [] communities;
        return 0;
    }

    uint32_t CopyPartition(CommunityPartition* destPartition, const CommunityPartition* sourcePartition) {
        destPartition->m_NodeLabels = new uint32_t[sourcePartition->m_NumNodes];
        if (!destPartition->m_NodeLabels) {
            printf("Error while allocating node labels.\n");
            return 1;
        }

        destPartition->m_CommunityIndices = new uint32_t[sourcePartition->m_NumCommunities];
        if (!destPartition->m_CommunityIndices) {
            printf("Error while allocating community indices.\n");
            return 1;
        }

        destPartition->m_Communities = new uint32_t[sourcePartition->m_NumCommunities + sourcePartition->m_NumNodes];
        if (!destPartition->m_Communities) {
            printf("Error while allocating communities.\n");
            return 1;
        }

        destPartition->m_InternalEdges = new uint32_t[sourcePartition->m_NumCommunities];
        if (!destPartition->m_InternalEdges) {
            printf("Error while allocating InternalEdges.\n");
            return 1;
        }

        destPartition->m_ExternalEdges = new uint32_t[sourcePartition->m_NumCommunities];
        if (!destPartition->m_ExternalEdges) {
            printf("Error while allocating InternalEdges.\n");
            return 1;
        }

        destPartition->m_NodeWCC = new double64_t[sourcePartition->m_NumNodes];
        if (!destPartition->m_NodeWCC) {
            printf("Error while allocating WCCs.\n");
            return 1;
        }

        memcpy(destPartition->m_NodeLabels, sourcePartition->m_NodeLabels, sizeof (uint32_t) * sourcePartition->m_NumNodes);
        memcpy(destPartition->m_CommunityIndices, sourcePartition->m_CommunityIndices, sizeof (uint32_t)*(sourcePartition->m_NumCommunities));
        memcpy(destPartition->m_Communities, sourcePartition->m_Communities, sizeof (uint32_t)*(sourcePartition->m_NumCommunities + sourcePartition->m_NumNodes));
        memcpy(destPartition->m_InternalEdges, sourcePartition->m_InternalEdges, sizeof (uint32_t)*(sourcePartition->m_NumCommunities));
        memcpy(destPartition->m_ExternalEdges, sourcePartition->m_ExternalEdges, sizeof (uint32_t)*(sourcePartition->m_NumCommunities));
        memcpy(destPartition->m_NodeWCC, sourcePartition->m_NodeWCC, sizeof (double64_t)*(sourcePartition->m_NumNodes));
        destPartition->m_NumNodes = sourcePartition->m_NumNodes;
        destPartition->m_NumCommunities = sourcePartition->m_NumCommunities;
        destPartition->m_WCC = sourcePartition->m_WCC;
        return 0;
    }


    uint32_t PrintPartition(const CGraph* graph, const CommunityPartition* partition, const char_t* fileName) {

        std::ofstream outFile;
        outFile.open(fileName);

        for (uint32_t i = 0; i < partition->m_NumCommunities; i++) {
            if (partition->m_CommunityIndices[i] != SCD_INVALID_COMMUNITY) {
                uint32_t* community = &partition->m_Communities[partition->m_CommunityIndices[i] + 1];
                uint32_t communitySize = partition->m_Communities[partition->m_CommunityIndices[i]];
                for (uint32_t j = 0; j < communitySize - 1; j++) {
                    outFile << graph->ReMap(community[j]) << " ";
                }
                outFile << graph->ReMap(community[communitySize - 1]) << std::endl;
            }
        }
        outFile.close();
        return 0;
    }


    void FreeResources(CommunityPartition* partition) {
        if (partition->m_NodeLabels != NULL) {
            delete [] partition->m_NodeLabels;
            partition->m_NodeLabels = NULL;
        }

        if (partition->m_CommunityIndices != NULL) {
            delete [] partition->m_CommunityIndices;
            partition->m_CommunityIndices = NULL;
        }

        if (partition->m_Communities != NULL) {
            delete [] partition->m_Communities;
            partition->m_Communities = NULL;
        }

        if (partition->m_InternalEdges != NULL) {
            delete [] partition->m_InternalEdges;
            partition->m_InternalEdges = NULL;
        }

        if (partition->m_ExternalEdges != NULL) {
            delete [] partition->m_ExternalEdges;
            partition->m_ExternalEdges = NULL;
        }

        if (partition->m_NodeWCC != NULL) {
            delete [] partition->m_NodeWCC;
            partition->m_NodeWCC = NULL;
        }
        partition->m_NumCommunities = 0;
        partition->m_NumNodes = 0;
        partition->m_WCC = 0.0;
    }


    uint32_t ImproveCommunities(const CGraph* graph, CommunityPartition* partition, uint32_t numThreads, uint32_t lookahead, const double64_t alfa ) {
        num_threads = numThreads;
        omp_set_num_threads(num_threads);
        printf("Maximum number of threads: %d\n", omp_get_max_threads());
        printf("Starting improvement from a partition with WCC: %f\n", partition->m_WCC / graph->GetNumNodes());
        CommunityPartition bestPartition;
        CopyPartition(&bestPartition, partition);

        uint32_t remainingTries = lookahead;
        bool improve = true;
        while(improve) {
            while (improve) {
                printf("\n");
                uint64_t initTime = StartClock();
                improve = false;
                printf("Starting improvement iteration ...\n");
                if (PerformImprovementStep(graph, partition, alfa)) {
                    printf("Error while performing an improvement step.\n");
                    return 1;
                }

                printf("New WCC: %f\n", partition->m_WCC / graph->GetNumNodes());
                printf("Best WCC: %f\n", bestPartition.m_WCC / graph->GetNumNodes());
                printf("Memory required by this iteration: %lu bytes \n", MeasureMemoryConsumption(partition) + MeasureMemoryConsumption(&bestPartition));

                if (partition->m_WCC - bestPartition.m_WCC > 0.0f) {
             //       if (((partition->m_WCC - bestPartition.m_WCC) / bestPartition.m_WCC) > 0.01f) {
                        remainingTries = lookahead;
              //      }
                    FreeResources(&bestPartition);
                    CopyPartition(&bestPartition, partition);
                } 


                printf("Iteration time: %lu ms\n", StopClock(initTime));
                if(remainingTries > 0) {
                    remainingTries--;
                    improve = true;
                }
            }
            printf("Trying to merge communities\n");
            MergeCommunities(graph, partition, alfa);
            printf("Merge: New WCC: %f\n", partition->m_WCC / graph->GetNumNodes());
            printf("Merge: Best WCC: %f\n", bestPartition.m_WCC / graph->GetNumNodes());
            if( partition->m_WCC - bestPartition.m_WCC > 0.0f ) {
                printf("Merging communities improved the partition\n");
            //    if (((partition->m_WCC - bestPartition.m_WCC) / bestPartition.m_WCC) > 0.01f) {
                    remainingTries = lookahead;
             //   }
                FreeResources(&bestPartition);
                CopyPartition(&bestPartition, partition);
                improve = true;
            }
        }

        FreeResources(partition);
        CopyPartition(partition, &bestPartition);
        FreeResources(&bestPartition);
        return 0;
    }
}
