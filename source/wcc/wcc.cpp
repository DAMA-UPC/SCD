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


        //            // v1.1
        //            while(list1 != endList1){
        //                while(list2 != endList2 && *list1 > *list2){
        //                    list2++;
        //                }
        //                
        //                if (*list1 == *list2){
        //                    triangles++;
        //                }                
        //                list1++;
        //            }


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

        // v 1.0
        //		uint32_t i = 0;
        //		uint32_t j = 0;
        //		uint32_t triangles = 0;
        //		while( i < size1 && j < size2 ) {
        //			if( list1[i] < list2[j] ) {
        //				i++;
        //			} 
        //			else if( list1[i] > list2[j] ) {
        //				j++;
        //			} else {
        //				triangles++;
        //				i++;
        //				j++;
        //			}
        //		}
        return triangles;
    }

    double64_t ComputeWCC(const CGraph * graph, const uint32_t * communities, const uint32_t numCommunities, const uint32_t* labelsIndices, const uint32_t * communitiesInvIndex, double64_t* wccs) {
        double64_t globalWCC = 0;
        #pragma omp parallel for schedule(static, 16) reduction(+:globalWCC)
        for (uint32_t i = 0; i < graph->GetNumNodes(); i++) {
            uint32_t communitySize = communitiesInvIndex[labelsIndices[communities[i]]];
            wccs[i] = ComputeWCC(graph, i, communities[i], communities, communitySize);
            globalWCC += wccs[i];
        }
        return globalWCC;
    }

    //TODO: Repassar per posar el codi d'interseccio de triangles

    double64_t ComputeWCC(const CGraph * graph, uint32_t node, uint32_t communityLabel, const uint32_t * communities, uint32_t communitySize) {
        uint32_t internalTriangles = 0;
        uint32_t internalTriangleDegree = 0;
        uint32_t triangleDegree = 0;
        uint32_t nodeId1 = node;
        const uint32_t* adjacencies1 = graph->GetNeighbors(nodeId1);
        uint32_t degree1 = graph->GetDegree(nodeId1);
        for (uint32_t k = 0; k < degree1; k++) {
            uint32_t nodeId2 = adjacencies1[k];
            uint32_t degree2 = graph->GetDegree(nodeId2);
            bool internal = false;
            if (communities[nodeId2] == communityLabel) {
                internal = true;
            }
            uint32_t i = 0;
            uint32_t j = 0;
            bool internalTriangleFound = false;
            bool triangleFound = false;
            const uint32_t* adjacencies2 = graph->GetNeighbors(nodeId2);
            while (i < degree1 && j < degree2) {
                if (adjacencies1[i] == adjacencies2[j]) {
                    uint32_t sharedNeighbor = adjacencies1[i];
                    if (internal && communities[sharedNeighbor] == communityLabel) {
                        internalTriangleFound = true;
                        internalTriangles++;
                    }
                    triangleFound = true;
                    i++;
                    j++;
                } else {
                    if (adjacencies1[i] < adjacencies2[j]) {
                        i++;
                    } else {
                        j++;
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
        if (graph->GetTotalTriangles(node) > 0 && communitySize > 1) {
            return ((internalTriangles / (double64_t) graph->GetTotalTriangles(node)) *
                    (triangleDegree / (double64_t) (triangleDegree - internalTriangleDegree + communitySize - 1)));
        }
        return 0.0;
    }
}


