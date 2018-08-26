#include "general_math.cuh"

bool isNan(double3 v)
{
	return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z);
}

char bool2int(bool b){ return b ? 1 : 0; }

int3 get_Nslc(double3 X, double3 dX){ return make_int3(int(X.x / dX.x), int(X.y / dX.y), int(X.z / dX.z)); }
int ind3D_to_ind(int3 I, int3 Sz)
{
	return I.x + Sz.x * (I.y + I.z * Sz.y);
}

bool is0(double3 v){ return (v.x == 0) && (v.y == 0) && (v.z == 0); }

double myRnd(void){ return rand()/double(RAND_MAX); }
double myRnd(double a, double b){
	if(a>b) swap(a,b);
	return myRnd()*(b-a)+a;
}
double3 vecByAngles(double phi, double tht )
{
	double ct = cos(tht);
	return make_double3(ct * cos(phi), ct * sin(phi), sin(tht));
}
double3 rndVec(double V)
{
	return vecByAngles(myRnd(-pi, pi), myRnd(-pi_d2, pi_d2))*V;
}
// this gauss is checked - it's really gauss
double gaussFnc(double x, double sgm, double x0)
{
	double b = (x-x0)/sgm;
	return exp(-b*b/2) / (sqrt(2*M_PI)*sgm);
}
double gaussRand(double sgm, double x0, double rng)
{
	rng *= sgm;
	double x;
	double y0 = 1.0/(sqrt(2*M_PI)*sgm); //y0 = gaussFnc(x0, sgm, x0); // max value
	double xl = x0-rng, xr = x0+rng;

	do{
		x = myRnd(xl, xr);
	}while(myRnd(0, y0) > gaussFnc(x, sgm, x0));

	return x;
}

bool almostEq(double a, double b, double _eps)
{
	return (b == 0 ? std::abs(a) : std::abs(a/b-1)) < _eps;
}
double epsDlt(double a, double b, double _eps)
{
	return b == 0 ? (a == 0 ? 0 : 1/_eps) : (std::abs(a) > std::abs(b) ? (a/b-1) : (b/a-1));
}
double float_part(double x)
{
    return x - int(x);
}

double cos_sin(double x)
{
    return sqrtf(1-x*x);
}

double sqr(double x)
{
	return x*x;
}

int sgn(double x, double _eps)
{
    if(x >= _eps){
    	return 1;
    } else if(x <= -_eps){
        return -1;
    } else {
        return 0;
    }
}
