#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <omp.h>
#include <algorithm>
#include <random>
#include <numeric>

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <unistd.h>
#include <new>
#include <iomanip>
#include <cuda_profiler_api.h>
#include <utility>

#include <math_helper.h>

using namespace std;

// --------------------------- math ----------------------------------

#define pi (3.14159265359)
#define pi_d2 (1.57079631679)
#define pi_m2 (6.28318530718)
//#define e (2.71828182846)
#define SQRT_2 (1.41421356237)
#define SQRT_pi (1.77245385090)
#define SQRT_2pi (2.50662827463)

// --------------------------- format ----------------------------------

#define PrintPres (6)
#define LOG_FILENAME ("comp.log")
#define TIME_UNIT (1000000) // MGz

#define PRINT_MODE_ID (200)
#define RAW_DATA_MODE (PRINT_MODE_ID)
#define FRAMES_DATA_MODE (PRINT_MODE_ID + 1)

// --------------------------- phys constants ---------------------------
#define SYS_EPS (1E-10)

#define RayTimeIND (100)
#define BaseRayType (RayTimeIND + 0)
#define PRayType (RayTimeIND + 1)
#define SRayType (RayTimeIND + 2)

/*
// -------------------------- element types -----------------------------

#define LINE_TYPE 1
#define TRIANGLE_TYPE 2
#define CUBE_TYPE 101
*/
// --------------------------- errors handling ---------------------------
#define SAY_IT (1)
#define ERROR_MSG (2)

#define FILE_ERROR_ID (10000)
#define CANT_OPEN_FILE_FOR_READING (FILE_ERROR_ID + 1)
#define CANT_OPEN_FILE_FOR_WRITING (FILE_ERROR_ID + 2)

#define WRONG_SURF_FILE_ID (FILE_ERROR_ID + 1000)
#define TOO_MANY_OVERLAPED_VERTICES (WRONG_SURF_FILE_ID + 1)
#define WRONG_VERTEX_INDEX (WRONG_SURF_FILE_ID + 2)

#define COMPUTATION_ERROR_ID (20000)
#define WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS (COMPUTATION_ERROR_ID + 1)
#define NO_REFLECTED_RAY_FOUND (COMPUTATION_ERROR_ID + 2)
#define WRONG_RAY_TYPE (COMPUTATION_ERROR_ID + 3)
#define LESS_2_MATERIALS (COMPUTATION_ERROR_ID + 4)

#define JUST_AN_ERROR_ID (30000)
#define NO_MODEL_NAME_SPECIFIED (JUST_AN_ERROR_ID + 1)

string error_handl_string, global_logFname = LOG_FILENAME;
