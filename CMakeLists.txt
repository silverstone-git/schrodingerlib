cmake_minimum_required(VERSION 3.16)
project(BallRolling LANGUAGES CXX)

set(CMAKE_AUTOMOC ON) # Required for Q_OBJECT and signals/slots

find_package(Qt6 REQUIRED COMPONENTS Widgets)

add_executable(BallRolling 
    main.cpp 
    BallCanvas.cpp 
    BallCanvas.h
)

target_link_libraries(BallRolling PRIVATE Qt6::Widgets)



# sudo apt update
# sudo apt install build-essential cmake qt6-base-dev
# mkdir build
# cd build
# cmake ..
# cmake --build .
# ./BallRolling
