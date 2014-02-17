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

#ifndef CGRAPH_H
#define CGRAPH_H

#include <assert.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <common/types.h>
#include <map>
#include <vector>
#include <set>

namespace scd {

		/**	@brief This class represents a graph.*/
class CGraph {
public:	
	CGraph();
	~CGraph();

	/** @brief Reads a graph from a file. The file must contain a list of
	*	undirected edges (an edge cannot appear twice). 
	*	The identifiers of the nodes have to bee between
	*	0 and N-1. The edges must be sorted by the first identifier first,
	*	and then the second.
	*	@param fileName The name of the file.
	*	@return 0 if the load was successful. 1 if there were errors.*/
	uint32_t 		Load(const char_t * fileName, uint32_t numThreads);

	/**	@brief Gets the number of nodes in the graph.
	*	@return The number of nodes.*/
	inline uint32_t 	GetNumNodes() const {
		return m_NumNodes;
	}

	/**	@brief Gets the number of edges in the graph.
	*	@brief The number of edges.*/
	inline uint32_t 	GetNumEdges() const {
			return m_NumEdges;
	}

	/**	@brief Gets the degree of a node.
	*	@param nodeId The identifier of the node.
	*	@return	The degree of the node.*/
	inline uint32_t 	GetDegree(uint32_t nodeId) const {
		assert(nodeId<m_NumNodes);
		return m_Nodes[nodeId].m_Degree;
	}

	/**	@brief Gets the neighbors of a node.
	*	@param nodeId The identifier of the node.
	*	@return	The neighbors of the node.*/
	inline const uint32_t * GetNeighbors(uint32_t nodeId) const {
		assert(nodeId<m_NumNodes);
		return &m_Adjacencies[m_Nodes[nodeId].m_AdjacencyIndex];
	}

	/** @brief Gets the identifier corresponding to the given node.
	*	@param nodeId The identifier of the node.
	*	@return The original identifier.*/
	inline uint32_t 	ReMap(uint32_t nodeId) const {
		assert(nodeId<m_NumNodes);
		return m_Map[nodeId];
	}

	/**	@brief Gets the total triangles of a node.
	*	@param nodeId The identifier of the node.
	*	@return	The total triangles of the node.*/
	inline uint32_t		GetTotalTriangles( uint32_t nodeId ) const {
		assert(nodeId<m_NumNodes);
		return m_TotalTriangles[nodeId];
	}

	/**	@brief Gets the clustering coefficient of the graph.
	 * 	@return The clustering coefficient.*/
	inline double64_t 	GetCC() const {
		return m_CC;
	}

	/** @brief Removes the edges that does not close any triangle, 
		and computes the number of triangles, the triangle degree and the clustering coefficient
		of each vertex.
		@param numThreads The number of threads to use.
		@return True if the removing was successful, false otherwise.*/
	uint32_t		RemoveEdgesNoTriangles(uint32_t numThreads);

private:
	uint32_t 		m_NumNodes; 		/**< @brief The number of nodes in the graph.*/
	uint32_t 		m_NumEdges;			/**< @brief The number of edges in the graph.*/
	Node* 			m_Nodes;			/**< @brief The array of nodes of the graph.*/
	uint32_t* 		m_Adjacencies;		/**< @brief The array of adjacencies of the graph.*/
	uint32_t* 		m_Map;				/**< @brief The map of internal identifiers to original identifiers.*/
	uint32_t* 		m_TotalTriangles;	/**< @brief The number of triangles a node belongs to.*/
	double64_t		m_CC;				/**< @brief The clustering coefficient of the graph.*/
};

}
#endif
