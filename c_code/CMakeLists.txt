cmake_minimum_required(VERSION 3.24)
project(c_code C)

set(CMAKE_C_STANDARD 99)

add_executable(c_code main.c step_1.c step_1.h step_2.c step_2.h)
