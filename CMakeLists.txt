# cmake_minimum_required(VERSION 3.16)
# project(BallRolling LANGUAGES CXX)
#
# set(CMAKE_AUTOMOC ON) # Required for Q_OBJECT and signals/slots
#
# find_package(Qt6 REQUIRED COMPONENTS Widgets)
#
# add_executable(BallRolling 
#     main.cpp 
#     BallCanvas.cpp 
#     BallCanvas.h
# )
#
# target_link_libraries(BallRolling PRIVATE Qt6::Widgets)
#
#
#
# # sudo apt update
# # sudo apt install build-essential cmake qt6-base-dev
# # mkdir build
# # cd build
# # cmake ..
# # cmake --build .
# # ./BallRolling

# ----------------------- New Non-Qt6 Project

cmake_minimum_required(VERSION 3.16)

# Set the C++ standard BEFORE the project or executable if possible
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Schrodingerlib LANGUAGES CXX)

# Find the package (VCPKG handles the heavy lifting)
find_package(Matplot++ CONFIG REQUIRED)

# Define the executable
# Syntax: add_executable(TargetName SourceFiles...)
add_executable(schrodinger_app main.cpp)

# Link the library
target_link_libraries(schrodinger_app PRIVATE Matplot++::matplot)
