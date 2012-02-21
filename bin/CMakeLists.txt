link_directories(${Boost_LIBRARY_DIRS})
include_directories(${CUDATHRUST_INCLUDE_DIR} ${CUDA_INCLUDE_DIR} ${TBB_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${CLUSTER_SOURCE_DIR}/src/)

cuda_add_executable(gpustat gpustat.cu)
target_link_libraries(gpustat ${CUDA_LIBRARY})

add_executable(graphstat graphstat.cpp)
target_link_libraries(graphstat ${Boost_LIBRARIES})

cuda_add_executable(clu clu.cpp)
target_link_libraries(clu ${CUDA_LIBRARY} ${TBB_LIBRARY} ${Boost_LIBRARIES})
target_link_libraries(clu cluster)

cuda_add_executable(dimacs dimacs.cpp)
target_link_libraries(dimacs ${CUDA_LIBRARY} ${TBB_LIBRARY} ${Boost_LIBRARIES})
target_link_libraries(dimacs cluster)

find_package(SDL)

if(SDL_FOUND)
	link_directories(${SDL_LIBRARY_DIR})
	include_directories(${SDL_INCLUDE_DIR})

	cuda_add_executable(visclu visclu.cpp)
	target_link_libraries(visclu ${CUDA_LIBRARY} ${TBB_LIBRARY} ${Boost_LIBRARIES} ${SDL_LIBRARY})
	target_link_libraries(visclu cluster)
endif()
