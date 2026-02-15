#pragma once

#include "const.cuh"

__host__ __device__ double cos_sin(double x);
__host__ __device__ double sqr(double x);
__host__ __device__ int sgn(double x, double _eps = SYS_EPS);
__host__ __device__ bool almostEq(double x, double y = 0, double _eps = SYS_EPS);
__host__ __device__ double epsDlt(double a, double b, double _eps = SYS_EPS);
__host__ __device__ double float_part(double x);

__host__ __device__ int3 get_Nslc(double3 X, double3 dX);
__host__ __device__ int ind3D_to_ind(int3 I, int3 Sz);
__host__ __device__ bool is0(double3 v);
__host__ __device__ char bool2int(bool b);
bool isNan(double3 v);

double myRnd(void);
double myRnd(double a, double b);
__host__ __device__ double3 vecByAngles(double phi, double tht = 0);
double3 rndVec(double V = 1);
__host__ __device__ double gaussFnc(double x, double sgm = 1, double x0 = 0);
double gaussRand(double sgm = 1, double x0 = 0, double rng = 5);
