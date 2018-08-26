#pragma once

#include "const.cuh"

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy);
int rightV(vector<double3>& v, double3 n, double3 v0);
