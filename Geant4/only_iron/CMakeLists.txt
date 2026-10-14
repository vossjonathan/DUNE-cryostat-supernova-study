#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)

project(B1 LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate build/compile_commands.json automatically
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


#----------------------------------------------------------------------------
# Geant4
#
find_package(Geant4 REQUIRED ui_all vis_all)


#----------------------------------------------------------------------------
# MARLEY / HepMC3 locations
#
set(MARLEY_DIR "/home/jonat/software/marley")

set(MARLEY_INCLUDE_DIR
    "${MARLEY_DIR}/include"
)

set(MARLEY_LIBRARY
    "${MARLEY_DIR}/build/lib/libMARLEY.so"
)

set(HEPMC3_LIBRARY
    "${MARLEY_DIR}/build/lib/libHepMC3.so"
)


# Find the directory which contains HepMC3/GenEvent.h
set(HEPMC3_INCLUDE_DIR
    "/home/jonat/software/marley/include/builtin"
)



#----------------------------------------------------------------------------
# Locate project source/header files
#
file(GLOB sources
    "${PROJECT_SOURCE_DIR}/src/*.cc"
)

file(GLOB headers
    "${PROJECT_SOURCE_DIR}/include/*.hh"
)


#----------------------------------------------------------------------------
# Create executable
#
add_executable(
    exampleB1
    exampleB1.cc
    ${sources}
    ${headers}
)


#----------------------------------------------------------------------------
# Include directories
#
target_include_directories(
    exampleB1
    PRIVATE

    # FDSteel headers
    "${PROJECT_SOURCE_DIR}/include"

    # Geant4 headers
    ${Geant4_INCLUDE_DIRS}

    # MARLEY headers
    "${MARLEY_INCLUDE_DIR}"

    # Needed by some MARLEY headers
    "${MARLEY_INCLUDE_DIR}/builtin"

    # HepMC3 headers
    "${HEPMC3_INCLUDE_DIR}"
)


#----------------------------------------------------------------------------
# Link libraries
#
target_link_libraries(
    exampleB1
    PRIVATE

    ${Geant4_LIBRARIES}

    "${MARLEY_LIBRARY}"
    "${HEPMC3_LIBRARY}"
)


# Allow executable to find MARLEY/HepMC3 shared libraries at runtime
set_target_properties(
    exampleB1
    PROPERTIES
    BUILD_RPATH "${MARLEY_DIR}/build/lib"
)


#----------------------------------------------------------------------------
# Copy macro/configuration scripts into build directory
#
set(EXAMPLEB1_SCRIPTS
    exampleB1.in
    exampleB1.out
    init_vis.mac
    run1.mac
    run2.mac
    vis.mac
)

foreach(_script ${EXAMPLEB1_SCRIPTS})
    configure_file(
        "${PROJECT_SOURCE_DIR}/${_script}"
        "${PROJECT_BINARY_DIR}/${_script}"
        COPYONLY
    )
endforeach()