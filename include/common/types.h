
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

#ifndef TYPES_H
#define TYPES_H

namespace scd {

	typedef bool bool_t;
	typedef unsigned char uchar_t;
	typedef char char_t;
	typedef short int uint16_t;
	typedef unsigned int uint32_t;
	typedef int int32_t;
	typedef long unsigned uint64_t;
	typedef float float32_t;
	typedef double double64_t;

	/** @brief This struct represents a node in the graph.*/
	struct Node {
		uint32_t m_Degree; 		/**< @brief The degree of the node.*/
		uint32_t m_AdjacencyIndex;	/**< @brief The index into the adjacency vector where the adjacencies lay*/
	};

	struct Edge {
		uint32_t head;
		uint32_t tail;

		bool operator<(const Edge & e) const
		{
			if(tail < e.tail)
			{
				return true;
			}
			else
			{
				if(tail == e.tail)
				{
					if(head < e.head)
					{
						return true;
					}
				}
			}
			return false;

		}
	};
};

#endif
