# 1. Minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

# 2. Project name and version
project(Varun VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT package
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# 4. Generate compile_commands.json for your Neovim LSP (clangd)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 5. Create the executable (Name of exe, then list of .cpp files)
add_executable(runme 
    main.cpp 
    polynomial/Polynomial.cpp
    integration/Leapfrog.cpp
    integration/Environment.cpp
    integration/Rocket.cpp
    heat/HeatEquation.cpp
    audio/AudioProcessor.cpp
)

target_include_directories(runme PRIVATE 
    .
    polynomial
    integration
    heat
    audio
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(runme PUBLIC ${ROOT_LIBRARIES})
