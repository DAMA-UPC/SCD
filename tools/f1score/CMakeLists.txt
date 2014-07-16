#SCD is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#SCD is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -std=c++11 -g3 -fno-rtti")	
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -std=c++11")	

AUX_SOURCE_DIRECTORY(./source SOURCE_FILES)
ADD_EXECUTABLE(f1score ${SOURCE_FILES})	
