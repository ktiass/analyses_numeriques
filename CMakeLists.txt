cmake_minimum_required(VERSION 3.25)
project(resolution_equation C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(resolution_equation
        control_fonction.c
        control_fonction.h
        equation_differentielle.c
        equation_differentielle.h
        equation_lineaire.c
        equation_lineaire.h
        interpolation_lineaire.c
        interpolation_lineaire.h
        main.c
        save_me.c
        save_me.h)
