
#include <fstream>
#include <common/types.h>
#include <graph/graph.h>
#include <iostream>
#include <math.h>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <queue>
#include <string.h>
#include <list>

using namespace scd;

std::map<uint32_t, uint32_t> id_map;

struct Subgraph
{
	unsigned int m_Size;
	std::set<unsigned int>* m_Adjacencies;

	Subgraph() {
		m_Size = 0;
        m_Adjacencies =  NULL;
	}

	~Subgraph() {
        Clear();
	}
    void Clear() {
        m_Size = 0;
        if(m_Adjacencies!=NULL) {
            delete [] m_Adjacencies;
            m_Adjacencies = NULL;
        }
    }
};


struct Statistics {
    unsigned int m_InternalTriangles;
    unsigned int m_TotalTriangles;
    unsigned int m_InternalDegree;
	unsigned int m_Diameter;	    // el diametre de la comunitat.
	unsigned int m_Bridges;			// nombre d'arestes pont dins la comunitat.
	unsigned int m_Size;		    // la mida de la comunitat.
    unsigned int m_ConnectedComponents;
    unsigned int m_TriangularComponents;
	double m_Wcc;
	double m_TriangleRatio;		    // trianglesIn / trianglesTotals.
	double m_TriangleDensity;		// trianglesIn / possibles trianglesIn.
    double m_EdgeDensity; 	        // per cada node, percentatge de nodes de la comunitat que son veins seus.
	double m_Conductance;			// conductança de la comunitat.
	double m_Tpr;				    // triangle participation ratio.
	double m_BridgeRatio;		    // nombre d'arestes pont dins la comunitat.
    double m_CC;
    double m_MinimumOverlapp;
    double m_MaximumOverlapp;

    Statistics() {
        Clear();
    }

    void Clear() {
        m_InternalTriangles = 0;
        m_TotalTriangles = 0;
        m_InternalDegree = 0;
        m_Wcc = 0.0;
        m_TriangleRatio = 0.0;		
        m_TriangleDensity = 0.0;		
        m_TriangularComponents = 0;
        m_ConnectedComponents = 0;
        m_EdgeDensity = 0.0; 	
        m_Conductance = 0.0;			
        m_Tpr = 0.0;				
        m_Diameter = 0;			
        m_Bridges = 0;			
        m_BridgeRatio = 0.0;		
        m_CC = 0.0;
        m_Size = 0;				
        m_MinimumOverlapp = 0.0;
        m_MaximumOverlapp = 0.0;
    }
};

struct Community {
    std::set<unsigned int> m_Nodes;
    Statistics m_Statistics;
    Subgraph m_Subgraph;
    int     m_Id;
    
    void Clear() {
        m_Nodes.clear();
        m_Statistics.Clear();
        m_Subgraph.Clear();
        m_Id = -1;
    }
};

std::vector<Community> communities;

bool ParseCommunity(std::ifstream& file, const scd::CGraph& graph, Community& community, int id) {
    community.Clear();
    std::string line;
	if(std::getline(file,line)) {
        community.m_Id = id;
        std::string delimiter(" ");
		size_t nextToken = 0;
		std::istringstream stream(line);
		unsigned int node;
		while( stream >> node ) {
			community.m_Nodes.insert(id_map.find(node)->second);
		}
        return true;
	}
    return false;
}

void ComputeSubgraph( const scd::CGraph& graph, Community& community ) {
	community.m_Subgraph.m_Size = community.m_Nodes.size();
	community.m_Subgraph.m_Adjacencies = new std::set<unsigned int>[community.m_Nodes.size()];
	std::map<unsigned int, unsigned int> mapa;
	unsigned int k = 0;
	for(std::set<unsigned int>::iterator it = community.m_Nodes.begin();it!=community.m_Nodes.end();it++,k++) {
		mapa.insert(std::pair<unsigned int, unsigned int>(*it,k));
	}
	k = 0;	
	for(std::set<unsigned int>::iterator it = community.m_Nodes.begin();it!=community.m_Nodes.end();it++,k++) {
		unsigned int node = *it;
		unsigned int degree = graph.GetDegree(node);
        const uint32_t * adjacencies = graph.GetNeighbors(node);
		for(unsigned int j = 0; j < degree;j++) {
			uint32_t neighbor = adjacencies[j];
			if(community.m_Nodes.find(neighbor) != community.m_Nodes.end()) {
				std::map<unsigned int,unsigned int>::iterator it2 = mapa.find(neighbor);
				if(it2==mapa.end()) {
					std::cout << "ERROR contructing subgraph " << std::endl;
				}
				community.m_Subgraph.m_Adjacencies[k].insert((*it2).second);
			}
		}
	}
}

void ComputeConnectedComponents(const CGraph& graph, Community& community) {
    bool * visited = new bool[community.m_Subgraph.m_Size];
    for(unsigned int j = 0; j < community.m_Subgraph.m_Size;++j) {
        visited[j] = false;
    }
    unsigned int connectedComponents = 0;
    for(unsigned int j = 0; j < community.m_Subgraph.m_Size;++j) {
        std::list<unsigned int> bfsList;
        if(!visited[j]) {
            connectedComponents++;
            visited[j] = true;
            bfsList.push_back(j);
            while(!bfsList.empty()) {
                unsigned int nextNode = bfsList.front();
                bfsList.pop_front();
                for(std::set<unsigned int>::iterator it = community.m_Subgraph.m_Adjacencies[nextNode].begin();it!=community.m_Subgraph.m_Adjacencies[nextNode].end();++it) {
                    unsigned int nextNeighbor = *it;
                    if(!visited[nextNeighbor]) {
                        visited[nextNeighbor] = true;
                        bfsList.push_back(nextNeighbor);
                    }
                }
            }
        }

    }
    delete [] visited;
    community.m_Statistics.m_ConnectedComponents = connectedComponents;
}



void ComputeTriangles( const CGraph& graph, Community& community, 
                       const unsigned int nodeId1,
                       const unsigned int nodeId2,
                       unsigned int& totalTriangles,
                       unsigned int& internalTriangles) {
	bool trianglePossible = true;
	if(community.m_Nodes.find(nodeId1)==community.m_Nodes.end() || community.m_Nodes.find(nodeId2)==community.m_Nodes.end()) {
		trianglePossible=false;
	}	
    totalTriangles = 0;
    internalTriangles = 0;
	unsigned int i = 0;
	unsigned int j = 0; 
    const uint32_t * adjacencies1 = graph.GetNeighbors(nodeId1);	
    const scd::uint32_t * adjacencies2 = graph.GetNeighbors(nodeId2);	
	while(i < graph.GetDegree(nodeId1) && j < graph.GetDegree(nodeId2)) {
			uint32_t node1 = adjacencies1[i]; 
			uint32_t node2 = adjacencies2[j]; 
			if( node1 == node2  ) {
				totalTriangles++;
				i++;
				j++;
				if(community.m_Nodes.find(node1)!=community.m_Nodes.end() && trianglePossible) {
                    internalTriangles++;
				}
			}
			else {
				if(adjacencies1[i] < adjacencies2[j]) {
					i++;
				} 
				else {
					j++;
				}
			}
	}
}

void ComputeTriangles( const CGraph& graph, Community& community) {
    std::set<unsigned int>::iterator it = community.m_Nodes.begin();
	for(;it!=community.m_Nodes.end();++it) {
		unsigned int node = *it;
		unsigned int degree=graph.GetDegree(node);
        const uint32_t * adjacencies = graph.GetNeighbors(node);
		for(unsigned int j = 0; j < degree;++j) {
			unsigned int totalTrianglesAux;
			unsigned int internalTrianglesAux;
			unsigned int neighbor = adjacencies[j];
			ComputeTriangles(graph,community,node,neighbor,totalTrianglesAux,internalTrianglesAux);
			community.m_Statistics.m_TotalTriangles+=totalTrianglesAux;
			community.m_Statistics.m_InternalTriangles+=internalTrianglesAux;
		}
	}
}


void ComputeTriangleDensity( const CGraph& graph, Community& community ) {
    ComputeTriangles(graph, community);
	community.m_Statistics.m_TriangleDensity = (double)community.m_Statistics.m_InternalTriangles/(double)(community.m_Nodes.size()*(community.m_Nodes.size()-1)*(community.m_Nodes.size()-2));
}

void ComputeTriangleRatio( const CGraph& graph, Community& community ) {
	community.m_Statistics.m_TriangleRatio = 0;
	if(community.m_Statistics.m_TotalTriangles>0){	
		community.m_Statistics.m_TriangleRatio = (double)community.m_Statistics.m_InternalTriangles/(double)community.m_Statistics.m_TotalTriangles;
	}	
}

void ComputeCC( const CGraph& graph, Community& community ) {
    std::set<unsigned int>::iterator it = community.m_Nodes.begin();
    unsigned int wEdges = 0;;
	for(;it!=community.m_Nodes.end();++it) {
		unsigned int node = *it;
		unsigned int degree=graph.GetDegree(node);
        unsigned int internalDegree = 0;
        const uint32_t * adjacencies = graph.GetNeighbors(node);
		for(unsigned int j = 0; j < degree;++j) {
			unsigned int neighbor = adjacencies[j];
			if(community.m_Nodes.find(neighbor)!=community.m_Nodes.end())  ++internalDegree;
		}
        wEdges += internalDegree*(internalDegree-1);
	}
    community.m_Statistics.m_CC = community.m_Statistics.m_InternalTriangles / (double)(wEdges);
}


void ComputeEdgeDensity( const CGraph& graph, Community& community ) {
    std::set<unsigned int>::iterator it = community.m_Nodes.begin();
	for(;it!=community.m_Nodes.end();++it) {
		unsigned int node = *it;
		unsigned int degree=graph.GetDegree(node);
        const uint32_t * adjacencies = graph.GetNeighbors(node);
		for(unsigned int j = 0; j < degree;++j) {
			unsigned int neighbor = adjacencies[j];
			if(community.m_Nodes.find(neighbor)!=community.m_Nodes.end()) community.m_Statistics.m_InternalDegree++;
		}
	}
    community.m_Statistics.m_EdgeDensity = community.m_Statistics.m_InternalDegree / (double)(community.m_Nodes.size()*(community.m_Nodes.size() - 1));
}

void ComputeBridges( const CGraph& graph, Community& community ) {
    ComputeSubgraph(graph,community);
	std::vector<bool> visited(community.m_Subgraph.m_Size,false);
	std::list<unsigned int> dfsList;
	std::list<unsigned int> stack;
	std::list<bool> booleanStack;
	unsigned int noParent = 0xffffffff;
	unsigned int numEdges = 0;
	for(unsigned int j = 0; j < community.m_Subgraph.m_Size;j++) {
		if( !visited[j] ) {
			unsigned int counter = 0;
			std::vector<unsigned int> labels(community.m_Subgraph.m_Size);
			std::vector<unsigned int> lowestLabels(community.m_Subgraph.m_Size);
			std::vector<unsigned int> parents(community.m_Subgraph.m_Size,noParent);
			dfsList.push_back(j);
			while(!dfsList.empty()){
				unsigned int nextNode = dfsList.back();
				dfsList.pop_back();
				if( !visited[nextNode] ) {
					stack.push_back(nextNode);
					booleanStack.push_back(true);
					visited[nextNode] = true;
					labels[nextNode] = counter;
					lowestLabels[nextNode] = counter;
					counter++;
					for(std::set<unsigned int>::iterator it = community.m_Subgraph.m_Adjacencies[nextNode].begin();
						it!=community.m_Subgraph.m_Adjacencies[nextNode].end();
						++it) {
						dfsList.push_back(*it);
						if( !visited[*it] ){
							parents[*it] = nextNode;
						}
					}
				} else if( parents[nextNode] != 0xffffffff && nextNode != parents[parents[nextNode]]) {
					stack.push_back(nextNode);
					booleanStack.push_back(false);
				}
			}
			while(!stack.empty()) {
				unsigned int nextNode = stack.back();		
				stack.pop_back();
				bool action = booleanStack.back();
				booleanStack.pop_back();
				if( action ) {
					if( parents[nextNode] != 0xffffffff ) {
						lowestLabels[parents[nextNode]] = std::min(lowestLabels[parents[nextNode]], lowestLabels[nextNode]);
						if(lowestLabels[nextNode] > labels[parents[nextNode]]) {
							community.m_Statistics.m_Bridges++;
						}
					}
				} else {
					for(std::set<unsigned int>::iterator it = community.m_Subgraph.m_Adjacencies[nextNode].begin();
						it!=community.m_Subgraph.m_Adjacencies[nextNode].end();
						++it) {
						if( parents[nextNode] != 0xffffffff && parents[nextNode] != *it) {
							lowestLabels[nextNode] = std::min(lowestLabels[nextNode],labels[*it]);
						}
					}
				}
			}
		}
		numEdges+=community.m_Subgraph.m_Adjacencies[j].size();
	}
}


void ComputeBridgeRatio( const CGraph& graph, Community& community ) {
    community.m_Statistics.m_BridgeRatio = (community.m_Statistics.m_Bridges*2) / (double)community.m_Statistics.m_InternalDegree;
}

void ComputeTPR( const CGraph& graph, Community& community ) {
		unsigned int numNodesWithTriangle = 0;
		for(std::set<unsigned int>::iterator it = community.m_Nodes.begin();it!=community.m_Nodes.end();it++) {
			unsigned int node = *it;
			unsigned int degree = graph.GetDegree(node);
			unsigned int outDegree = 0;
			unsigned int numTriangles = 0;
            const uint32_t * adjacencies = graph.GetNeighbors(node);
			for(unsigned int j = 0; j < degree; j++) {
				if(community.m_Nodes.find(adjacencies[j]) != community.m_Nodes.end()) {
					unsigned int totalTriangles;
					unsigned int internalTriangles;
					ComputeTriangles(graph,community,node, adjacencies[j], totalTriangles, internalTriangles );
					if( internalTriangles > 0 ) {
						numTriangles++;
					}
				}
			}
			if( numTriangles > 0 ) {
				numNodesWithTriangle++;
			}
		}
		community.m_Statistics.m_Tpr = numNodesWithTriangle/ (double)community.m_Nodes.size();
}

static	unsigned int SelectMax( const std::vector<int>& distances ) {
		unsigned int size = distances.size();
		assert(size > 0);
		int max = distances[0];
		unsigned int maxNode = 0;
		for( unsigned int i = 1; i < size; ++i ) {
			if( distances[i] > max ) {
				maxNode = i;
				max = distances[i];
			}
		}
		return maxNode;
	}

void ComputeDiameter(const CGraph& graph, Community& community) {
	std::vector<int> distances(community.m_Subgraph.m_Size, 0xffffffff);
	int maxExcentricity = 0;
	unsigned int numIterations = sqrt(community.m_Subgraph.m_Size);	
	for( unsigned int k = 0; k < numIterations; k++ ) {
		std::vector<bool> visited(community.m_Subgraph.m_Size, false);
		std::list<unsigned int> bfsList;
		unsigned int node;
		if( k == 0 ) { node = rand() % community.m_Subgraph.m_Size; } 
		else { node = SelectMax(distances); } 
		bfsList.push_back(node);
		std::list<int> levelList;
		levelList.push_back(-1);
		while(!bfsList.empty()) {
			unsigned int nextNode = bfsList.front();
			bfsList.pop_front();
			int level = levelList.front() + 1;
			levelList.pop_front();
			distances[nextNode] = std::min( distances[nextNode], level );
			maxExcentricity = std::max(maxExcentricity, level);
			for(std::set<unsigned int>::iterator it = community.m_Subgraph.m_Adjacencies[nextNode].begin();it!=community.m_Subgraph.m_Adjacencies[nextNode].end();it++) {
				unsigned int nextNeighbor = *it;
				if(!visited[nextNeighbor]) {
					visited[nextNeighbor] = true;
					bfsList.push_back(nextNeighbor);
					levelList.push_back(level);
				}
			}
		}
	}
	community.m_Statistics.m_Diameter = static_cast<unsigned int>(maxExcentricity);
}

void ComputeConductance(const CGraph& graph, Community& community) {
	unsigned int outDegree = 0;
	unsigned int totalDegree = 0;
	std::set<unsigned int>::iterator it = community.m_Nodes.begin();
	for(;it!=community.m_Nodes.end();it++) {
		unsigned int node = *it;
		unsigned int degree = graph.GetDegree(node);
        const uint32_t * adjacencies = graph.GetNeighbors(node);
		totalDegree += degree;
		for(unsigned int j = 0; j < degree;j++) {
			unsigned int neighbor = adjacencies[j];
			if(community.m_Nodes.find(neighbor) == community.m_Nodes.end()) {
				outDegree++;
			}
		}
	}
	community.m_Statistics.m_Conductance = 0;
	if(totalDegree>0) {
		community.m_Statistics.m_Conductance = outDegree/(double)totalDegree;
	}
}

void RemoveNoTriangles( const CGraph& graph, Community& community, Community& outCommunity ) {
    unsigned int numNodesWithTriangle = 0;
    for(std::set<unsigned int>::iterator it = community.m_Nodes.begin();it!=community.m_Nodes.end();it++) {
        unsigned int node = *it;
        unsigned int degree = graph.GetDegree(node);
        unsigned int outDegree = 0;
        unsigned int numTriangles = 0;
        const uint32_t * adjacencies = graph.GetNeighbors(node);
        for(unsigned int j = 0; j < degree; j++) {
            if(community.m_Nodes.find(adjacencies[j]) != community.m_Nodes.end()) {
                unsigned int totalTriangles;
                unsigned int internalTriangles;
                ComputeTriangles(graph,community,node, adjacencies[j], totalTriangles, internalTriangles );
                if( internalTriangles > 0 ) {
                    numTriangles++;
                }
            }
        }
        if( numTriangles > 0 ) {
            outCommunity.m_Nodes.insert(node);
        }
    }
}

void ComputeTriangularComponents( const CGraph& graph, Community& community) {
    Community auxCommunity;
    RemoveNoTriangles(graph,community,auxCommunity);
    ComputeSubgraph(graph,auxCommunity);
    ComputeConnectedComponents(graph,auxCommunity);
    community.m_Statistics.m_TriangularComponents = auxCommunity.m_Statistics.m_ConnectedComponents;
}


void ComputeOverlapp( const CGraph& graph, Community& community, std::map<int, std::vector< int > >& nodeCommunities ) {
    std::set<int> alreadyChecked;
    int minOverlapp = community.m_Nodes.size();
    int maxOverlapp = 0;
    for( std::set<unsigned int>::iterator it = community.m_Nodes.begin(); it != community.m_Nodes.end(); ++it ) {
        const std::vector< int >& coms = nodeCommunities[*it];
        for( std::vector< int >::const_iterator it2 = coms.begin(); it2 != coms.end(); ++it2 ) {
            if( (*it2) != community.m_Id && alreadyChecked.find(*it2) == alreadyChecked.end() ) {
                std::set<unsigned int> aux;
                const std::set< unsigned int >& other = communities[*it2].m_Nodes;
                //std::cout << other.size() << std::endl;
    //            std::cout << community.m_Nodes.size() << std::endl;
                std::set_intersection( other.begin(), other.end(), community.m_Nodes.begin(), community.m_Nodes.end(), std::inserter(aux,aux.begin()));
                minOverlapp = minOverlapp > aux.size() ? aux.size() : minOverlapp;
                maxOverlapp = maxOverlapp < aux.size() ? aux.size() : maxOverlapp;
                alreadyChecked.insert( *it2 );
            } 
        } 
    }
    community.m_Statistics.m_MinimumOverlapp = minOverlapp / (double)community.m_Nodes.size();
    community.m_Statistics.m_MaximumOverlapp = maxOverlapp / (double)community.m_Nodes.size();
}

void ComputeStatistics( const CGraph& graph, Community& community ) {
    ComputeTriangleDensity(graph,community);
    ComputeTriangleRatio(graph,community);
    ComputeEdgeDensity(graph,community);
    ComputeCC(graph,community);
    ComputeBridges(graph,community);
    ComputeBridgeRatio(graph,community);
    ComputeTPR(graph,community);
    ComputeConductance(graph,community);
    ComputeDiameter(graph,community);
    ComputeTriangularComponents(graph,community);
}

#define CHECK_ARGUMENT_STRING(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = argv[index+1]; \
			} else { \
				std::cout << "Invalid options" << std::endl; \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_FLOAT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atof(argv[index+1]); \
			} else { \
				std::cout << "Invalid options" << std::endl; \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_INT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atoi(argv[index+1]); \
			} else { \
				std::cout << "Invalid options" << std::endl; \
				return 1;\
			}\
		}

#define CHECK_FLAG(index, option,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
		}

void PrintUsage() {
    std::cout << "communityAnalyzer -f <network_file> -p <partition_file> -o <output_file>" << std::endl;
}

void PrintCommunity( std::ofstream& file, const Community& community ) {
    file << community.m_Id << "\t"
         << community.m_Statistics.m_TriangleDensity << "\t"
         << community.m_Statistics.m_CC << "\t"
         << community.m_Statistics.m_EdgeDensity << "\t"
         << community.m_Statistics.m_TriangleRatio << "\t"
         << community.m_Nodes.size() << "\t"
         << community.m_Statistics.m_Diameter << "\t"
         << community.m_Statistics.m_Tpr << "\t"
         << community.m_Statistics.m_BridgeRatio << "\t"
         << community.m_Statistics.m_Conductance << "\t"
         << community.m_Statistics.m_TriangularComponents << "\t"
         << community.m_Statistics.m_MinimumOverlapp << "\t"
         << community.m_Statistics.m_MaximumOverlapp << "\t"
         << std::endl;
}

void PrintFileHeader( std::ofstream& file) {
    file << "Id\t"
         << "TriangleDensity\t"
         << "CC\t"
         << "EdgeDensity\t"
         << "TriangleRatio\t"
         << "Size\t"
         << "Diameter\t"
         << "Tpr\t"
         << "BridgeRatio\t"
         << "Conductance\t"
         << "TriangularComponents\t"
         << "MinimumOverlapp\t"
         << "MaximumOverlapp\t"
         << std::endl;
}

int main(int argc, char ** argv) {
    char* graphFileName = NULL;
    char* partitionFileName = NULL;
    char* outputFileName = NULL;
    bool graphFileNameSet = false;
    bool partitionFileNameSet = false;
    bool outputFileNameSet = false;

	for( unsigned int i = 1; i < argc; i++) {
		CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
		CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
		CHECK_ARGUMENT_STRING(i, "-p", partitionFileName, partitionFileNameSet)
	}

    if( !graphFileNameSet || !outputFileNameSet || !partitionFileNameSet ) {
        PrintUsage();
        exit(1);
    }

    CGraph graph;
	if(graph.Load(graphFileName,1)!=0) {
		std::cout << "ERROR: Unable to load graph" << std::endl;
	}

    const uint32_t* map_array = graph.GetMap();
    for( uint32_t i = 0; i < graph.GetNumNodes(); ++i) {
        id_map[map_array[i]] = i;
    }

    std::ofstream outputFile;
    outputFile.open(outputFileName);
    PrintFileHeader(outputFile);
    std::ifstream partitionFile;
    partitionFile.open(partitionFileName);
    std::map< int, std::vector<int> > nodeCommunityIndex;
    int numParsed = 0;
    Community community;
    while(ParseCommunity(partitionFile,graph,community, numParsed)) {
        communities.push_back(community);
        for( std::set< unsigned int >::iterator it = community.m_Nodes.begin(); it != community.m_Nodes.end(); ++it ) {
            nodeCommunityIndex[*it].push_back(community.m_Id);
        }
        ++numParsed;
        if( numParsed % 1000 == 0 ) {
            std::cout << "Parsed " << numParsed << " communities" << std::endl;
        }
    }

    int numComputed = 0;
    for( int i = 0; i < communities.size(); ++i ) {
        if(communities[i].m_Nodes.size() > 2) {
            ComputeStatistics(graph,communities[i]);
            ComputeOverlapp( graph, communities[i], nodeCommunityIndex );
            PrintCommunity(outputFile,communities[i]);
        }
        ++numComputed;
        if( numComputed % 1000 == 0 ) {
            std::cout << "Computed statistics of " << numComputed << " communities" << std::endl;
        }
    }
    partitionFile.close();
    outputFile.close();
	return 0;
}
