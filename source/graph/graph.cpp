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
#include <graph/graph.h>
#include <list>
#include <map>
#include <wcc/wcc.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <cstdlib>

namespace scd {

			/**	@brief Compares two unsigned integers.
			 * 	@param e1 Void pointer to the first unsigned integer.	
			 * 	@param e2 Void pointer to the second unsigned integer.	
			 *	@return -1 if e1 goes before e2. 1 if e1 goes after e2. 0 if e1 and e2 are equal.*/
	static int 	Compare_Ids(const void* e1, const void* e2) {
		uint32_t id1 = *(uint32_t*)e1;
		uint32_t id2 = *(uint32_t*)e2;
		if( id1 < id2 ) return -1;
		if( id2 < id1 ) return 1;
		return 0;
	}
        

        
        

	CGraph::CGraph() :
	m_NumNodes(0),
	m_NumEdges(0),
	m_Nodes(NULL),
	m_Adjacencies(NULL),
	m_Map(NULL),
	m_TotalTriangles(NULL)
	{

	}

	CGraph::~CGraph() {

		if( m_Nodes!=NULL ) {
			delete [] m_Nodes;
			m_Nodes = NULL;
		}

		if( m_Adjacencies!=NULL ) {
			delete [] m_Adjacencies;
			m_Adjacencies = NULL;
		}

		if( m_Map != NULL ) {
			delete []  m_Map;
			m_Map = NULL;
		}
		
		if( m_TotalTriangles != NULL ) {
			delete []  m_TotalTriangles;
			m_TotalTriangles = NULL;
		}
	}

	uint32_t CGraph::Load(const char_t * fileName, uint32_t numThreads) {

		printf(	"Graph: Loading Graph\n" );
		std::ifstream inFile;
		inFile.open((const char *)fileName);
		if(!inFile) {
			printf( "Graph: Error Openning Graph File\n" );
			return 1;
		}
             

		printf( "Graph: Relabeling nodes ...\n" );
		timeval time;
		gettimeofday(&time, NULL);
		uint64_t initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		std::map<uint32_t, uint32_t>* mapa = new std::map<uint32_t,uint32_t>();
		if( !mapa ) {
			printf( "\t Graph: Error allocating mapa\n" );
			return 1;
		}
		uint32_t index = 0;
		m_NumEdges = 0;
		uint32_t node1;
		while( inFile >> node1 ) {
			uint32_t node2;
			inFile >> node2;

			if(!mapa->count(node1)) {
				mapa->insert(std::pair<uint32_t,uint32_t>(node1,index));
				index++;
			}

			if(!mapa->count(node2)) {
				mapa->insert(std::pair<uint32_t,uint32_t>(node2,index));
				index++;
			}
			m_NumEdges++;
		}
		gettimeofday(&time, NULL);
		uint64_t endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		printf("Graph: Nodes relabeled in %lu ms\n", endTime - initTime);


		printf( "Graph: Reading degrees ...\n" );
		//We set the file cursor to the beginning.
		inFile.close();
		inFile.open((const char *)fileName);
		gettimeofday(&time, NULL);
		initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		//Allocate space for nodes and initialize de degree.
		m_NumNodes = index;
		m_Nodes = new Node[m_NumNodes];
		if( !m_Nodes ) {
			printf( "Graph: Error allocating nodes\n" );
			return 1;
		}
		for( uint32_t i = 0; i < m_NumNodes; i++ ) {
			m_Nodes[i].m_Degree = 0;
		}

		//Compute the degree of each node.
		while( inFile >> node1 ) {
			uint32_t node2;
			inFile >> node2;
			m_Nodes[(*mapa->find(node1)).second].m_Degree++;
			m_Nodes[(*mapa->find(node2)).second].m_Degree++;
		}

		//Computing the adjacency indices, average degree and maximum Degree.
		float32_t averageDegree = 0.0f;
		float32_t maxDegree = 0.0f;
		uint32_t currentAdjacencyIndex = 0;
		for( uint32_t i = 0; i < m_NumNodes; i++ ) {
			m_Nodes[i].m_AdjacencyIndex = currentAdjacencyIndex;
			currentAdjacencyIndex += m_Nodes[i].m_Degree;
			averageDegree += m_Nodes[i].m_Degree;
			if( m_Nodes[i].m_Degree > maxDegree ) {
				maxDegree = m_Nodes[i].m_Degree;
			}
		}
		averageDegree /= (m_NumNodes);
		gettimeofday(&time, NULL);
		endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		printf("\t Graph: Degrees read in %lu ms\n", endTime - initTime);

		printf( "Graph: Reading adjacencies ...\n" );
		gettimeofday(&time, NULL);
		initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		//We set the file cursor to the beginning.
		inFile.close();
		inFile.open((const char *)fileName);
		
		m_Adjacencies = new uint32_t[m_NumEdges*2];
		if( !m_Adjacencies ) {
			printf( "Graph: Error allocating adjacencies\n" );
			return 1;
		}
		uint32_t* counters = new uint32_t[m_NumNodes];
		for( uint32_t i = 0; i < m_NumNodes; i++ ) {
			counters[i] = 0;
		}

		//Filling adjacencies
		while( inFile >> node1 ) {
			uint32_t node2;
			inFile >> node2;
			uint32_t tail = (*mapa->find(node1)).second;
			uint32_t head = (*mapa->find(node2)).second;
			assert(counters[tail]<m_Nodes[tail].m_Degree);
			assert(counters[head]<m_Nodes[head].m_Degree);
			m_Adjacencies[m_Nodes[tail].m_AdjacencyIndex + counters[tail]] = head;
			m_Adjacencies[m_Nodes[head].m_AdjacencyIndex + counters[head]] = tail;
			counters[tail]++;
			counters[head]++;
		}
		delete [] counters;
		inFile.close();
		gettimeofday(&time, NULL);
		endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		printf("\t Graph: Adjacencies read in %lu ms\n", endTime - initTime);

		printf( "Graph: Sorting adjacencies ...\n" );
		gettimeofday(&time, NULL);
		initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
                
		//Sorting adjacencies.
                #pragma omp parallel for schedule(static, 16)
		for(uint32_t i = 0; i < m_NumNodes; i++ ) {
			qsort(&m_Adjacencies[m_Nodes[i].m_AdjacencyIndex], m_Nodes[i].m_Degree, sizeof(uint32_t),Compare_Ids);
		}
		gettimeofday(&time, NULL);
		endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		printf("\t Graph: Adjacencies sorted in %lu ms\n", endTime - initTime);
		
		printf( "Graph: Filling map array ...\n" );
		gettimeofday(&time, NULL);
		initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		//Filling the map array
		m_Map = new uint32_t[m_NumNodes];
		for(std::map<uint32_t,uint32_t>::iterator it = mapa->begin();it!=mapa->end();it++) {
			m_Map[(*it).second] = (*it).first;
		}	

		delete mapa;
		gettimeofday(&time, NULL);
		endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
		printf("\t Graph: Map array filled in %lu ms\n", endTime - initTime);

		printf( "Graph: Graph Loaded\n" );
		printf( "Graph: Number of Nodes: %u\n", m_NumNodes );
		printf( "Graph: Number of Edges: %u\n", m_NumEdges );
		//printf( "Graph: Clustering coefficient: %f\n", m_CC ); //Not available here
		printf( "Graph: Average Degree: %f\n", averageDegree );
		printf( "Graph: Maximum Degree: %f\n", maxDegree );
		printf( "Graph: Memory \n" );
		printf( "..............\n" );
		uint64_t memNodes = (char*)&m_Nodes[m_NumNodes] - (char*)&m_Nodes[0];
		uint64_t memEdges = (char*)&m_Adjacencies[m_NumEdges*2] - (char*)&m_Adjacencies[0];
		uint64_t memMap   = (char*)&m_Map[m_NumNodes-1] - (char*)&m_Map[0];
		uint64_t memTotalTriangles =  (char*)&m_TotalTriangles[m_NumNodes-1] - (char*)&m_TotalTriangles[0];
		printf( "%-16s %-10lu Bytes\n", "Nodes:", memNodes );
		printf( "%-16s %-10lu Bytes\n", "Adjacencies:", memEdges );
		printf( "%-16s %-10lu Bytes\n", "Map:", memMap );
		printf( "%-16s %-10lu Bytes\n", "TotalTriangles:", memTotalTriangles );
		printf( "%-16s %-10lu Bytes\n", "Total:", memNodes + memEdges + memMap + memTotalTriangles );
		printf( "..............\n" );
		return 0;
	}


        static int compareInt (const void * a, const void * b)        
        {
            if ( *(int*)a >  *(int*)b ) return 1;
            if ( *(int*)a <  *(int*)b ) return -1;
            if ( *(int*)a == *(int*)b ) return 0;   
        }
        
        
	uint32_t	CGraph::RemoveEdgesNoTriangles( uint32_t numThreads) {
		omp_set_num_threads(numThreads);
		m_TotalTriangles = new uint32_t[m_NumNodes];
		if( !m_TotalTriangles ) {
			printf("Error allocating total triangles\n");
			return 1;
		}
		m_CC = 0;
		uint32_t numEdgesRemoved = 0;
		uint32_t newAdjacencyIndex = 0;
		uint32_t* edgesTriangles  = new uint32_t[m_NumEdges*2];
                
		#pragma omp parallel for schedule(dynamic, 32)
		for(uint32_t i = 0; i < m_NumNodes; i++ ) {
                        uint32_t edgesTrianglesIndex = m_Nodes[i].m_AdjacencyIndex;
			m_TotalTriangles[i] = 0;			
			uint32_t* adjacencyList1 = &m_Adjacencies[m_Nodes[i].m_AdjacencyIndex];
			for(uint32_t j = 0; j < m_Nodes[i].m_Degree; j++) {
				uint32_t* adjacencyList2 = &m_Adjacencies[m_Nodes[adjacencyList1[j]].m_AdjacencyIndex];
				if( i < adjacencyList1[j]) { 
					uint32_t triangles =  Intersect(adjacencyList1, m_Nodes[i].m_Degree, adjacencyList2, m_Nodes[adjacencyList1[j]].m_Degree);
					edgesTriangles[edgesTrianglesIndex] = triangles;
                                        
                                        uint32_t k = (uint32_t*) bsearch(&i, adjacencyList2, 
                                                       m_Nodes[adjacencyList1[j]].m_Degree, sizeof(uint32_t), compareInt)
                                                     - adjacencyList2;                                        
                                        assert(k !=  m_Nodes[adjacencyList1[j]].m_Degree); // "ERROR when computing triangles."                                        
                                        edgesTriangles[m_Nodes[adjacencyList1[j]].m_AdjacencyIndex + k ] = triangles;
				}
				edgesTrianglesIndex++;
			}
		}

		uint32_t edgesTrianglesIndex = 0;
		for(uint32_t i = 0; i < m_NumNodes; i++ ) {
			m_TotalTriangles[i] = 0;
			uint32_t  newDegree = 0;
			uint32_t* tempAdjacencies = new uint32_t[m_Nodes[i].m_Degree];
			uint32_t* adjacencyList1  = &m_Adjacencies[m_Nodes[i].m_AdjacencyIndex];
			for(uint32_t j = 0; j < m_Nodes[i].m_Degree; j++) {
				uint32_t triangles = edgesTriangles[edgesTrianglesIndex++];
                                if( triangles > 0 ) {
                                    tempAdjacencies[newDegree] = adjacencyList1[j];
                                    newDegree++;
                                } else {
                                    numEdgesRemoved++;
                                }
                                m_TotalTriangles[i] += triangles;
			}
                        
			memcpy(&m_Adjacencies[newAdjacencyIndex], tempAdjacencies, sizeof(uint32_t)*newDegree);
			m_Nodes[i].m_Degree = newDegree;
			m_Nodes[i].m_AdjacencyIndex = newAdjacencyIndex;
			newAdjacencyIndex+=newDegree;
			delete [] tempAdjacencies;
			uint32_t auxPossibleTriangles = m_Nodes[i].m_Degree*(m_Nodes[i].m_Degree - 1);
			if( auxPossibleTriangles > 0 ) {
				m_CC += m_TotalTriangles[i] / (double64_t)auxPossibleTriangles;
			}                        
		}
                
                delete [] edgesTriangles;
		m_CC /= m_NumNodes;
                std::cout << "m_cc inicial: " << m_CC << std::endl;
		m_NumEdges-=(numEdgesRemoved/2);
		return 0;
	}
}
