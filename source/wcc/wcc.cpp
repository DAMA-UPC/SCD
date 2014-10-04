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

#include <wcc/wcc.h>
#include <omp.h>
#include <communities/communities.h>

namespace scd {

    
    uint32_t Intersect(/*const*/ uint32_t* list1, const uint32_t size1, /*const*/ uint32_t* list2, const uint32_t size2) {
        uint32_t triangles = 0;

        uint32_t* endList1 = list1 + size1;
        uint32_t* endList2 = list2 + size2;


        // v2.0
        while (list1 != endList1 && list2 != endList2) {
            if (*list1 < *list2) {
                list1++;                
            } else if (*list1 > *list2) {
                list2++;                
            } else { //(*list1 == *list2){ //triangle found
                triangles++;
                list1++;
                list2++;
            }
        }

        return triangles;
    }
    
    double64_t ComputeWCC(const CGraph * graph, const double64_t alfa, const uint32_t * communities, 
            const uint32_t numCommunities, const uint32_t* labelsIndices, 
            const uint32_t * communitiesInvIndex, double64_t* wccs) {    
        double64_t globalWCC = 0;
        
        #pragma omp parallel for schedule(static, 8) reduction(+:globalWCC)
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            uint32_t communitySize = communitiesInvIndex[labelsIndices[communities[i]]];
            wccs[i] = ComputeWCC(graph, alfa, i, communities[i], communities, communitySize);
            globalWCC += wccs[i];
        }
        return globalWCC;
    }

    
    double64_t ComputeWCC(const CGraph * graph, const double64_t alfa, uint32_t node, uint32_t communityLabel, 
            const uint32_t * communities, uint32_t communitySize) {
        
        uint32_t        internalTriangles      = 0;
        uint32_t        internalTriangleDegree = 0;
        uint32_t        triangleDegree         = 0;
        uint32_t        node1                  = node;
        const uint32_t* adjacencies1           = graph->GetNeighbors(node1);
        uint32_t        degree1                = graph->GetDegree(node1);
        
        if (communitySize <= 2 ||graph->GetTotalTriangles(node) == 0) {
            return 0.0;
        }
        

        //while(adjacencies1 < endList1){
        for (uint32_t k = 0; k < degree1; k++) {
            uint32_t        nodeId2  = adjacencies1[k];
            uint32_t        degree2  = graph->GetDegree(nodeId2);                        
            bool            internal = (communities[nodeId2] == communityLabel);
            bool            internalTriangleFound = false;
            bool            triangleFound         = false;
            const uint32_t* adjacencies2          = graph->GetNeighbors(nodeId2);
            
            uint32_t*  currentNode1     = (uint32_t*) adjacencies1;
            uint32_t*  currentNode2     = (uint32_t*) adjacencies2;
            uint32_t*  endAdjacencies1  = (uint32_t*) adjacencies1 + degree1;
            uint32_t*  endAdjacencies2  = (uint32_t*) adjacencies2 + degree2;
            
            while (currentNode1 != endAdjacencies1 && currentNode2 != endAdjacencies2){
                if (*currentNode1 == *currentNode2){
                    uint32_t sharedNeighbor = *currentNode1;
                    if (internal && communities[sharedNeighbor] == communityLabel) {
                        internalTriangleFound = true;
                        internalTriangles++;
                    }
                    triangleFound = true;
                    currentNode1++;
                    currentNode2++;
                }else  if(*currentNode1 < *currentNode2){
                    while(*currentNode1 < *currentNode2 && currentNode1 < endAdjacencies1){
                        currentNode1++;
                    }
                }else{
                    while(*currentNode1 > *currentNode2 && currentNode2 < endAdjacencies2){
                        currentNode2++;
                    }
                }
            }
            
            if (internalTriangleFound) {
                internalTriangleDegree++;
            }
            if (triangleFound) {
                triangleDegree++;
            }
        }
        
        return ((internalTriangles / (double64_t) graph->GetTotalTriangles(node)) *
               (triangleDegree / (double64_t) (triangleDegree + alfa*(communitySize - 1 - internalTriangleDegree))));
    }
}


