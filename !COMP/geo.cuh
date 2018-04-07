/*
 * Space.cuh
 *
 *  Created on: Jul 7, 2016
 *      Author: polyachyadt
 */

//#ifndef GEO_CUH_
//#define GEO_CUH_

#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <fstream>
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cstring>
#include <sys/stat.h>
#include <omp.h>
#include <string>
#include <sstream>
#include <unistd.h>
#include <new>
#include <algorithm>
#include <iomanip>
#include <cuda_profiler_api.h>
#include <random>
#include <numeric>
#include <utility>

#include <math_helper.h>

using namespace std;

class Ray;
class Triangle;
class Surface;
class Params;

// --------------------------- phys constants ---------------------------
//#define Cp 800
//#define Cs 500
#define SYS_EPS (1E-10)
#define LOG_FILENAME "comp.log"
//#define Nray 10
//#define MapSize 10000
//#define Tmax (MapSize/Cs)
//#define Amin 0.00001

#define RayTimeIND 100
#define PRayType (RayTimeIND + 1)
#define SRayType (RayTimeIND + 2)

// --------------------------- errors handling ---------------------------
#define FILE_ERROR_ID (10000)
#define SAY_IT FILE_ERROR_ID
#define CANT_OPEN_FILE_FOR_READING (FILE_ERROR_ID + 1)
#define CANT_OPEN_FILE_FOR_WRITING (FILE_ERROR_ID + 2)

#define WRONG_MSH_FILE_ID (FILE_ERROR_ID + 1000)
#define TOO_BIG_NODE_IND (WRONG_MSH_FILE_ID + 1)
#define TOO_BIG_ELEMENT_IND (WRONG_MSH_FILE_ID + 2)

#define MATH_ERROR_ID (20000)
#define WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS (MATH_ERROR_ID + 1)
#define NO_REFLECTED_RAY_FOUND (MATH_ERROR_ID + 2)

template<typename T>
int CHECK(int n, T s, string logFname = LOG_FILENAME);

template<typename T>
void SAY_LOG(T s, string logFname = LOG_FILENAME);

// --------------------------- some math ---------------------------
double cos_sin(double x);
double sqr(double x);
int sgn(double x, double _eps = 0);
bool almostEq(double x, double y = 0, double _eps = SYS_EPS);

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy);
int rightV(vector<double3>& v, double3 n, double3 v0);

class Params{
public:
	double Cp, Cs;
	double eps;
	double Tmax, Amin;
	int Nrays;
	long long int n_alive_rays;
	long long int n_total_rays;

	string msh_filename;

	Params(void):
		Cp(0), Cs(0), eps(SYS_EPS), Tmax(0), Amin(0), Nrays(0), n_alive_rays(0), n_total_rays(0), msh_filename("")
	{}
	~Params(void){}

	int load_from_file(string filename);
	int print(ostream &output, string spr = "\n");
	int save_to_file(string filename);
};

Params prm;

class Ray{
public:
    int type;
    double c,A,t;
    double3 r,v;

    Ray():
            type(0), c(0), A(0), t(0)
    {
            v = r = make_double3(0,0,0);
    }
    Ray(const double3 _r, const double3 _ve, const int _type = PRayType, const double _A = 1, const double _t = 0):
            type(_type), A(_A), t(_t)
    {
    	this->r = _r;
        this->c = (this->type == PRayType ? prm.Cp : prm.Cs);
        this->v = _ve * this->c;
    }
     ~Ray(){}

    void print(ostream &output, string spr = "\n-----\n");
    void move(Surface& srf);
    pair<double, int> find_collision(Surface& srf);
};

class Triangle{
public:
    double3 r[3];
    double3 n;

    Triangle(){ this->r[0] = this->r[1] = this->r[2] = this->n = make_double3(0,0,0); }
    Triangle(const double3 r1, const double3 r2, const double3 r3){
            this->r[0] = r1;
            this->r[1] = r2;
            this->r[2] = r3;
            this->n = this->getNorm();
    }

    double3 getNorm(void) const;
    int sg(double3 rx) const;
    int isInside(double3 rx) const;
};

class  Surface{
public:
	Triangle *polygons;
	int size;

    Surface(int _n = 0);

    void resize_clear(int _n);
    int load_from_file(string filename);
};

class World{
	Params prms;
	Surface srf;
};

template<typename T>
void stp(T str);

template<typename T>
string vectorToStr(vector<T> v, string sp = ";");

template<typename T>
void printVector(ostream &output, vector<T> v, string sp1="(", string sp2=";", string sp3=")\n");

template<typename T>
T sumVector(vector<T> v);

template <typename T>
string toString(T val);

template<typename T>
T fromString(const string& s);

string toLower(string s);
string toUpper(string s);
void printD3(ostream &output, double3 r, string sp1="(", string sp2=";", string sp3=")\n");

