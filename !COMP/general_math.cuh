#pragma once

#include "const.cuh"

double cos_sin(double x);
double sqr(double x);
int sgn(double x, double _eps = SYS_EPS);
bool almostEq(double x, double y = 0, double _eps = SYS_EPS);
double epsDlt(double a, double b, double _eps = SYS_EPS);
double float_part(double x);

int3 get_Nslc(double3 X, double3 dX);
int ind3D_to_ind(int3 I, int3 Sz);
bool is0(double3 v);
char bool2int(bool b);
bool isNan(double3 v);

double myRnd(void);
double myRnd(double a, double b);
double3 vecByAngles(double phi, double tht = 0);
double3 rndVec(double V = 1);
double gaussFnc(double x, double sgm = 1, double x0 = 0);
double gaussRand(double sgm = 1, double x0 = 0, double rng = 5);
