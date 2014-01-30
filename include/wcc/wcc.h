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
	

#ifndef WCC_H
#define WCC_H


#include <common/types.h>
#include <graph/graph.h>
namespace scd 
{

			/**	@brief Computes the size of the intersection between two arrays.
			 * 	@param[in] list1 The first array.
			 * 	@param[in] size1 The size of the first array.
			 * 	@param[in] list2 The second array.	
			 * 	@param[in] size2 The size of the second array.
			 * 	@return The size of the intersection.*/
	uint32_t 	Intersect( const uint32_t* list1, const uint32_t size1, const uint32_t* list2, const uint32_t size2 );

			/**	@brief Computes the WCC of a node against a community.
			 * 	@param[in] graph The graph.
			 * 	@param[in] communities The assignment of nodes to communities.
			 * 	@param[in] labelsIndices The array of indexes of labels into the community inverse index.
			 * 	@param[in] communitiesInvIndex The community inverse index.
			 * 	@param[in] wccs An array where the WCCs of the nodes will be stored.
			 * 	@return The WCC of the node against the community.*/
	double64_t 	ComputeWCC(const CGraph * graph, const uint32_t * communities, const uint32_t numCommunities,  const uint32_t* labelsIndices, const uint32_t * communitiesInvIndex, double64_t* wccs);

			/**	@brief Computes the WCC of a node against a community.
			 * 	@param[in] graph The graph.
			 * 	@param[in] node The node.
			 * 	@param[in] communityLabel The label of the community to test against.
			 * 	@param[in] communities The assignment of nodes to communities.
			 * 	@param[in] numCommunities The number of communities.
			 * 	@param[in] labelsIndices The array of indexes of labels into the community inverse index.
			 * 	@param[in] communitiesInvIndex The community inverse index.
			 * 	@return The WCC of the node against the community.*/
	double64_t 	ComputeWCC(const CGraph * graph, uint32_t node, uint32_t communityLabel, const uint32_t * communities, uint32_t communitySize );
}

#endif



