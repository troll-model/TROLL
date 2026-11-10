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

# --- Test target ---
# check.sh $1 = nombre de pas de temps, défaut 20 dans le script
foreach(n 10 100 1000 10000)
    add_custom_target(sanity${n}
        COMMAND ${CMAKE_SOURCE_DIR}/sanity/check.sh ${n}
        COMMENT "Running sanity check — ${n} iterations (bit-by-bit comparison)"
        DEPENDS ${PROJECT_NAME}
    )
endforeach()

# Make 'test' depend on 'sanity' with 100 iterations (and any other tests in the future)
add_custom_target(test DEPENDS sanity100)

# --- Clean target ---
set_property(DIRECTORY APPEND PROPERTY
    ADDITIONAL_CLEAN_FILES
        "${CMAKE_BINARY_DIR}/base_results"
        "${CMAKE_BINARY_DIR}/test_runs"
)

