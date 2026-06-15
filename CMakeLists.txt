cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 14)

project(quasar)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 /WX-")
  set(CMAKE_CXX_FLAGS_DEBUG "/Zi /Od")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-math-errno -fno-trapping-math")
endif()

FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

# OpenMP — Apple Clang needs libomp from Homebrew; hint its location
if(APPLE)
  execute_process(COMMAND brew --prefix libomp
                  OUTPUT_VARIABLE HOMEBREW_LIBOMP
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
  if(HOMEBREW_LIBOMP)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP}/lib/libomp.dylib")
  endif()
endif()

find_package(OpenMP COMPONENTS CXX)

add_executable(quasar
  ${CMAKE_SOURCE_DIR}/src/Quasar.cpp
  ${CMAKE_SOURCE_DIR}/src/Data.cpp
  ${CMAKE_SOURCE_DIR}/src/Geno.cpp
  ${CMAKE_SOURCE_DIR}/src/Utils.cpp
  ${CMAKE_SOURCE_DIR}/src/Phi.cpp 
  ${CMAKE_SOURCE_DIR}/src/Residualise.cpp 
  ${CMAKE_SOURCE_DIR}/src/ModelFit.cpp 
  ${CMAKE_SOURCE_DIR}/src/ScoreTest.cpp 
  ${CMAKE_SOURCE_DIR}/src/QTLMappingUtils.cpp
)
target_include_directories(quasar PRIVATE ${CMAKE_SOURCE_DIR}/src)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found — building with parallel score tests")
  target_link_libraries(quasar OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found — building single-threaded")
endif()

set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")

target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/cxxopts/include/)
target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/eigen-3.4.0/)
target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/brent/)

enable_testing()

add_executable(quasar_tests EXCLUDE_FROM_ALL tests/test_smoke.cpp)
target_include_directories(quasar_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_compile_definitions(quasar_tests PRIVATE
    EXAMPLE_DATA_DIR="${CMAKE_SOURCE_DIR}/example"
    QUASAR_EXECUTABLE="$<TARGET_FILE:quasar>"
    BUILD_DIR="${CMAKE_BINARY_DIR}"
)

add_test(NAME smoke_lm COMMAND quasar_tests "LM smoke test")
add_test(NAME smoke_nb_glm COMMAND quasar_tests "NB-GLM smoke test")
add_test(NAME smoke_lmm COMMAND quasar_tests "LMM smoke test")
