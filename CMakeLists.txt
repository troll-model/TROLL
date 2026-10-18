cmake_minimum_required(VERSION 3.13)

project(TROLL VERSION 4.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}, build type ${CMAKE_BUILD_TYPE}")

# --- Options -----------------------------------------------------------------
option(TROLL_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (-Werror); ON in CI" OFF)
option(TROLL_EXTRA_WARNINGS "Add -Wextra -Wshadow -Wconversion (advisory: several hundred warnings today)" OFF)
set(TROLL_SANITIZE "" CACHE STRING "Semicolon-separated sanitizers, e.g. address;undefined")
set(TROLL_DEFINES "" CACHE STRING
    "Compile-time model options, e.g. \"Output_ABC;TROLL_NO_WATER\" (see include/troll_defines.hpp)")

# --- Dependencies --------------------------------------------------------------
find_package(GSL REQUIRED)
message(STATUS "GSL ${GSL_VERSION} found in ${GSL_INCLUDE_DIRS}")

# --- Executable ----------------------------------------------------------------
add_executable(${PROJECT_NAME}
    src/troll.cpp
    src/species.cpp
    src/crown.cpp
    src/tree.cpp
    src/lookup_tables.cpp
    src/input.cpp
    src/memory.cpp
    src/output.cpp
    src/simulation.cpp
    src/internal/lookup_tables_cache.cpp
    src/internal/param_registry.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/internal
)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas)

# Warnings are ON (the former global "-w" silently disabled all of them).
target_compile_options(${PROJECT_NAME} PRIVATE -Wall)
if(TROLL_EXTRA_WARNINGS)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wextra -Wshadow -Wconversion)
endif()
if(TROLL_WARNINGS_AS_ERRORS)
    target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
endif()

if(TROLL_DEFINES)
    message(STATUS "TROLL model options: ${TROLL_DEFINES}")
    target_compile_definitions(${PROJECT_NAME} PRIVATE ${TROLL_DEFINES})
endif()

if(TROLL_SANITIZE)
    string(REPLACE ";" "," _troll_sanitizers "${TROLL_SANITIZE}")
    message(STATUS "Sanitizers: ${_troll_sanitizers}")
    target_compile_options(${PROJECT_NAME} PRIVATE -fsanitize=${_troll_sanitizers} -fno-omit-frame-pointer -g)
    target_link_options(${PROJECT_NAME} PRIVATE -fsanitize=${_troll_sanitizers})
endif()

# --- Regression checks (see sanity/README.md) --------------------------------------
set(_troll_check ${CMAKE_SOURCE_DIR}/sanity/check.sh --troll $<TARGET_FILE:${PROJECT_NAME}>)

# 20 days from bare ground on the 4 ha example (fast smoke test)
add_custom_target(sanity
    COMMAND ${_troll_check} --case example --nbiter 20
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: example, 20 days (bit-exact)" USES_TERMINAL)

# 1000 days from bare ground (replaces the missing check_1000.sh)
add_custom_target(sanity1000
    COMMAND ${_troll_check} --case example --nbiter 1000
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: example, 1000 days (bit-exact)" USES_TERMINAL)

# One year on an established 1 ha stand: adult growth, competition, mortality, treefall
add_custom_target(sanity-stand
    COMMAND ${_troll_check} --case stand --nbiter 365
    DEPENDS ${PROJECT_NAME} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/sanity
    COMMENT "Sanity check: established stand, 1 year (bit-exact)" USES_TERMINAL)

add_custom_target(check DEPENDS sanity sanity-stand)
