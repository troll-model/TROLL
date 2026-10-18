cmake_minimum_required(VERSION 3.10.0)
project(TROLL_LINT LANGUAGES NONE)

# --- Linting targets ---
# Find cpplint and cppcheck
find_program(CPPLINT_PATH cpplint)
find_program(CPPCHECK_PATH cppcheck)

# cpplint target
if(CPPLINT_PATH)
    add_custom_target(lint-style
        COMMAND ${CPPLINT_PATH} ../include/*.hpp ../src/*.cpp
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running cpplint for code style checks"
    )
else()
    message(WARNING "cpplint not found, lint-style target will not be available")
endif()

# cppcheck target
if(CPPCHECK_PATH)
    add_custom_target(lint-static
        COMMAND ${CPPCHECK_PATH}
        ARGS --enable=all --inconclusive --std=c++17 --suppress=missingIncludeSystem ../src/
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running cppcheck for static analysis"
    )
    # Gating check used in CI: warning-level findings fail the build.
    # uninitMemberVar is suppressed for now (82 members of Tree/Species without
    # initializers); remove the suppression once they are initialized.
    add_custom_target(lint-static-gate
        COMMAND ${CPPCHECK_PATH}
        ARGS --enable=warning,portability --inline-suppr --std=c++17
             --suppress=uninitMemberVar --suppress=missingIncludeSystem
             -I../include -I../include/internal --error-exitcode=1 --quiet ../src/
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running cppcheck (gating: warnings are errors)"
    )
else()
    message(WARNING "cppcheck not found, lint-static target will not be available")
endif()

# Combined lint target
add_custom_target(lint
    DEPENDS lint-style lint-static
    COMMENT "Running all linting tools (style + static analysis)"
)