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


#ifndef COMMUNITIES_H
#define COMMUNITIES_H

#include <graph/graph.h>
#include <list>

namespace scd {

#define SCD_SINGLETON	0xffffffff

	/**	@brief This struct is a tuple formed by a node id and a clustering coefficient.*/
	struct NodeClustering {
		uint32_t m_NodeId;
		double64_t m_CC;
		uint32_t m_Degree;
	};

	struct CommunityPartition {
		uint32_t* 	m_NodeLabels;		/**< @brief The labels of the communities each node belongs to.*/
		uint32_t*	m_CommunityIndices;	/**< @brief The array of indices for each label into the community array.*/
		uint32_t*	m_Communities;		/**< @brief The communities.*/
		uint32_t*	m_InternalEdges;	/**< @brief The number of internal edges of each community.*/
		uint32_t*	m_ExternalEdges;	/**< @brief The number of external edges of each community.*/
		double64_t*	m_NodeWCC;			/**< @brief The WCC of the nodes.*/
		uint32_t	m_NumCommunities;	/**< @brief The number of communities.*/
		uint32_t	m_NumNodes;			/**< @brief The number of nodes.*/
		double64_t	m_WCC;				/**< @brief The WCC of this partition.*/
	};

			/**	@brief Initializes a partition structure with an initial partition. 
			 * 	@param[in] graph A pointer to the graph.	
			 * 	@param[out] partition The partition structure to initializes. 
			 *	@return 0 if the computation was successful. 1 if there were errors.*/
	uint32_t 	InitializeSimplePartition( const CGraph* graph, CommunityPartition* partition );
		
			/**	@brief Frees the resources used by the partition.
			 *	@param[out] partition The partition to free.*/
	void		FreeResources( CommunityPartition* partition );

			/**	@brief Prints the communities into a file.
			 * 	@param[in] graph The graph.
			 * 	@param[in] partition The partition to print.
			 * 	@param[in] fileName The name of the file where the communities will be print.
			 * 	@return 0 if the execution was successful. 1 otherwise.*/
	uint32_t	PrintPartition( const CGraph* graph, const CommunityPartition* partition, const char_t* fileName);


			/**	@brief Improves the quality of the communities.
			 * 	@param[in] graph The graph.
			 * 	@param[out] partition The partition to improve.
			 * 	@param[in] numThreads The number of threads.
			 * 	@return 0 if the execution was successful. 0 otherwise.*/ 
	uint32_t	ImproveCommunities( const CGraph* graph, CommunityPartition* partition, uint32_t numThreads );

			/**	@brief Copies a partition into another partition.
			 * 	@param[out] destPartition The destination partition.
			 * 	@param[in] sourcePartition The source partition.
			 * 	@resutl 0 if the copy was sucecssful. 1 otherwise.*/
	uint32_t	CopyPartition( CommunityPartition* destPartition, const CommunityPartition* sourcePartition);
}


#endif


