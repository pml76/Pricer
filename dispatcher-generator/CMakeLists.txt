cmake_minimum_required(VERSION 3.12)
project(dispgen)

set(SOURCES dispatcher-generator.cpp)

add_executable(dispgen ${SOURCES})