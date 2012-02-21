#Copyright 2011, Bas Fagginger Auer.
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
cmake_minimum_required(VERSION 2.6)
project(CLUSTER)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CLUSTER_SOURCE_DIR})

find_package(CUDA REQUIRED)
find_package(CUDAThrust REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost COMPONENTS iostreams program_options REQUIRED)

#set(CMAKE_CXX_FLAGS "-O2 -g -Wall -Wextra -Wshadow -ansi -pedantic -DTIME")
#set(CMAKE_CXX_FLAGS "-O2 -g -Wall -Wextra -Wno-long-long -Wno-unused-parameter -ansi -DTIME")
set(CMAKE_CXX_FLAGS "-O3 -Wall -ffast-math -DNDEBUG -DLEAN")
#set(CUDA_NVCC_FLAGS -ccbin gcc-4.4 -Xptxas -dlcm=cg)
#set(CUDA_NVCC_FLAGS -ccbin gcc-4.4)
#set(CUDA_NVCC_FLAGS -Xptxas -dlcm=cg)

add_subdirectory(${CLUSTER_SOURCE_DIR}/src/)
add_subdirectory(${CLUSTER_SOURCE_DIR}/bin/)

