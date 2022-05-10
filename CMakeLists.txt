cmake_minimum_required(VERSION 3.21)
project(maxrank CXX)

set(CMAKE_CXX_STANDARD)
set(SRC ./src)
set(LIB ./lib)
set(INC ./include)

set(RTREECORE_HEADERS
        ${INC}/collection.h
        ${INC}/rtree.h ${INC}/rstartree.h ${INC}/rnode.h ${INC}/rentry.h ${INC}/point.h ${INC}/hypercube.h ${INC}/line2d.h ${INC}/virtualRnode.h ${INC}/Query.h ${INC}/tgs.h
        ${INC}/mainmem.h
        ${INC}/filemem.h
        ${INC}/search.h ${INC}/psdraw.h
        ${INC}/iomeasure.h)
set(RTREECORE_SOURCES
        ${LIB}/collection.cc
        ${LIB}/rtree.cc ${LIB}/rstartree.cc ${LIB}/rnode.cc ${LIB}/rentry.cc ${LIB}/point.cc ${LIB}/hypercube.cc ${LIB}/line2d.cc ${LIB}/virtualRnode.cc ${LIB}/Query.cc ${LIB}/tgs.cc
        ${LIB}/mainmem.cc
        ${LIB}/filemem.cc
        ${LIB}/search.cc ${LIB}/psdraw.cc
        ${LIB}/iomeasure.cc)

set(MYQUADTREE_HEADERS
        ${INC}/QuadTree.h ${INC}/S3.h ${INC}/param.h)
set(MYQUADTREE_SOURCES
        ${SRC}/myQuadTree.cpp ${SRC}/QuadTree.cpp ${SRC}/S3.cc ${SRC}/param.cc)

add_library(rtreecore ${RTREECORE_SOURCES} ${RTREECORE_HEADERS})
target_include_directories(rtreecore PRIVATE ${INC})

add_executable(myQuadTree ${MYQUADTREE_SOURCES} ${MYQUADTREE_HEADERS})
target_include_directories(myQuadTree PRIVATE ${INC})
target_link_libraries(myQuadTree PUBLIC rtreecore)
