cmake_minimum_required(VERSION 3.10.0)

project(TROLL VERSION 4.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
MESSAGE(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}")

# --- Build target ---
# Find GSL package
find_package(GSL REQUIRED)
if(NOT ${GSL_FOUND})
    MESSAGE(FATAL_ERROR "The GSL library was not found")
endif()
MESSAGE(STATUS "GSL found")
MESSAGE(STATUS "   - version ${GSL_VERSION}")
MESSAGE(STATUS "   - included files located in ${GSL_INCLUDE_DIRS}")


add_compile_options(-w)
add_executable(${PROJECT_NAME}
    src/troll.cpp
    src/species.cpp
    src/field.cpp
    src/tree.cpp
    src/lookup_tables.cpp
    src/input.cpp
    src/memory.cpp
    src/output.cpp
    src/simulation.cpp
    src/internal/lookup_tables_cache.cpp
    src/internal/param_registry.cpp
)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # Public headers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/internal  # Internal headers
)

target_compile_options(${PROJECT_NAME} PRIVATE -O3 -Wall)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas)

# sanity check
add_custom_target(sanity
    COMMAND ${CMAKE_SOURCE_DIR}/sanity/check.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Running sanity check (bit-by-bit comparison)"
)

# sanity check — 1000 iterations
add_custom_target(sanity1000
    COMMAND ${CMAKE_SOURCE_DIR}/sanity/check_1000.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Running sanity check — 1000 iterations (bit-by-bit comparison)"
)

# Make 'test' depend on 'sanity' (and any other tests in the future)
add_custom_target(test DEPENDS sanity)


