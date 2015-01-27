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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)
PROJECT(SCD CXX)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -pg -fopenmp -DPROFILE ")	
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG")	
#SET(CMAKE_VERBOSE_MAKEFILE ON)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

ADD_SUBDIRECTORY(./tools/wcc)
ADD_SUBDIRECTORY(./tools/f1score)
ADD_SUBDIRECTORY(./tools/selector)
ADD_SUBDIRECTORY(./tools/cc)

INCLUDE_DIRECTORIES(./include)
FILE( GLOB_RECURSE SOURCE_FILES "source/*" )
ADD_EXECUTABLE(scd ${SOURCE_FILES})	
