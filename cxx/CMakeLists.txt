cmake_minimum_required(VERSION 3.10)

set(EXEC_NAME DNA_bind)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set the project name
project(${EXEC_NAME}  VERSION 0.1 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#compiler flags
if (MSVC)
	SET(CMAKE_CXX_FLAGS_RELEASE "/O2")
	SET(CMAKE_C_FLAGS_RELEASE "/O2")
elseif(APPLE) #assuming clang
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -march=native")
	set(CMAKE_CXX_FLAGS_DEBUG "-g3 -ggdb")
	set(CMAKE_CXX_FLAGS_RELEASE "-flto=full -O2 -Wall -Wextra -Wpedantic -Wno-unused-parameter -ftree-vectorize -march=native")
	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
else() #assuing gcc
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -march=native")
	set(CMAKE_CXX_FLAGS_DEBUG "-g3 -ggdb")
	set(CMAKE_CXX_FLAGS_RELEASE "-flto -O2 -Wall -Wextra -Wpedantic -Wno-unused-parameter -ftree-vectorize -msse2 -mfpmath=sse -march=native")
	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")	
endif()
       
# sources
set(SRCS ./src/main.cpp ./src/bindmodel.cpp ./src/derivedmodels.cpp ./src/dna.cpp)

# add the executable
add_executable(${EXEC_NAME} ${SRCS})

#headers
target_include_directories(${EXEC_NAME}  PRIVATE ./src/include/dnabind ./src/include/utilities ./src/include/data)
target_include_directories(${EXEC_NAME}  SYSTEM PRIVATE ./vendor/armadillo ./vendor) #SYSTEM keyword suppresses warnings from this files

if (MSVC)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT supported OUTPUT error)
	if( supported )
	    message(STATUS "IPO / LTO enabled")
	    set_property(TARGET ${EXEC_NAME}  PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
	else()
	    message(STATUS "IPO / LTO not supported: <${error}>")
	endif()
	option(CMAKE_USE_WIN32_THREADS_INIT "using WIN32 threads" ON)
elseif (APPLE)
	set_target_properties(${EXEC_NAME} PROPERTIES LINK_FLAGS_RELEASE -flto=full)
else()
	set_target_properties(${EXEC_NAME} PROPERTIES LINK_FLAGS_RELEASE -flto)
endif()

find_package(Threads REQUIRED)
target_link_libraries(${EXEC_NAME} PUBLIC Threads::Threads)

#Python - note you may need to set the version to your specific version for which Nupack4.0 has been installed
#find_package(Python 3.8.10 REQUIRED) #<---- modern version gives empty string for ${PYTHON_INCLUDE_DIRS}$, defaulting to old version below
set(PY_VER 3.8)
find_package(PythonInterp ${PY_VER} REQUIRED) #<------ proceeding with old version
find_package(PythonLibs ${PY_VER} REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
message(STATUS "found python headers at " ${PYTHON_INCLUDE_DIRS})
target_link_libraries(${EXEC_NAME} PRIVATE ${PYTHON_LIBRARIES})
