#pragma once

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
class Material;

// --------------------------- format ----------------------------------

#define PrintPres 6
#define LOG_FILENAME "comp.log"

// --------------------------- phys constants ---------------------------
//#define Cp 800
//#define Cs 500
#define SYS_EPS (1E-10)
//#define Nray 10
//#define MapSize 10000
//#define Tmax (MapSize/Cs)
//#define Amin 0.00001

#define RayTimeIND 100
#define BaseRayType (RayTimeIND + 0)
#define PRayType (RayTimeIND + 1)
#define SRayType (RayTimeIND + 2)

// -------------------------- element types -----------------------------

#define LINE_TYPE 1
#define TRIANGLE_TYPE 2
#define CUBE_TYPE 101

// --------------------------- errors handling ---------------------------
#define FILE_ERROR_ID (10000)
#define SAY_IT FILE_ERROR_ID
#define CANT_OPEN_FILE_FOR_READING (FILE_ERROR_ID + 1)
#define CANT_OPEN_FILE_FOR_WRITING (FILE_ERROR_ID + 2)

#define WRONG_SURF_FILE_ID (FILE_ERROR_ID + 1000)
#define TOO_MANY_OVERLAPED_VERTICES (WRONG_SURF_FILE_ID + 1)

#define COMPUTATION_ERROR_ID (20000)
#define WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS (COMPUTATION_ERROR_ID + 1)
#define NO_REFLECTED_RAY_FOUND (COMPUTATION_ERROR_ID + 2)
#define WRONG_RAY_TYPE (COMPUTATION_ERROR_ID + 3)
#define LESS_2_MATERIALS (COMPUTATION_ERROR_ID + 4)

#define JUST_AN_ERROR_ID (30000)
#define NO_MODEL_NAME_SPECIFIED (JUST_AN_ERROR_ID + 1)

template<typename T>
//int CHECK(int n, T s, string logFname = LOG_FILENAME);
int CHECK(int n, T s);

template<typename T>
void SAY_LOG(T s);

// --------------------------- some math ---------------------------
double cos_sin(double x);
double sqr(double x);
int sgn(double x, double _eps = SYS_EPS);
bool almostEq(double x, double y = 0, double _eps = SYS_EPS);

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy);
int rightV(vector<double3>& v, double3 n, double3 v0);
Material *getMatByInd(int Nmat, Material* mat, int ind0);

class Material{
public:
	double Cp, Cs;

	Material(double _cp = 0, double _cs = 0):
		Cp(_cp), Cs(_cs)
	{}
	bool isEq(const Material* m2, const double _eps = SYS_EPS) const;

	string toStr(string spr1 = "", string spr2 = "");
	void print(ostream &output, string spr1 = "", string spr2 = "\n");
};

class Params{
public:
	double eps;
	double Tmax, Amin;
	int Nrays;
	long long int n_alive_rays;
	long long int n_total_rays;

	string model_name, tet_filename, prm_filename, material_filename;

	vector<string> paramFHead;

	Params(void):
		eps(SYS_EPS), Tmax(0), Amin(0), Nrays(0), n_alive_rays(0), n_total_rays(0), model_name("sys_tst")
	{
		this->paramFHead = {"Nrays","Tmax","Amin","eps"};
	}
	~Params(void){}

	int load_from_file(string filename);
	void print(ostream &output, string spr = "\n");
	int save_to_file(string filename);
};

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
    Ray(const double3 _r, const double3 _v, const int _type = BaseRayType, const double _A = 1, const double _t = 0):
            type(_type), A(_A), t(_t)
    {
    	this->r = _r;
        this->v = _v;
        this->c = (this->type == BaseRayType ? 1 : length(_v));
    }
     ~Ray(){}

    void print(ostream &output, string spr = "\n-----\n");
    string toStr(const string spr = "\n-----\n");
    void move(Surface* srf, Params *prm);
    pair<double, int> find_collision(Surface* srf, Params* prm);
};

class Triangle{
public:
    double3 r[3];
    int i_vrtx[3]; // r[i] == nodes[i_vrtx[i]]
    double3 n; // n = [(r1-r0)x(r2-r0)]

    Material* mat[2];
    // mat[1] is in the side where normal-vector points, mat[0] is in the other (back) side

    Triangle();
    Triangle(const double3 r1, const double3 r2, const double3 r3);
    Triangle(const double3* r_new);
    Triangle(const double3 *vertex, const int* ind);
    Triangle(const double3 *vertex, const int i0, const int i1, const int i2);

    void clear_params();
    double3 getNorm(void) const;
    int sg(const double3 rx, const double _eps = SYS_EPS) const;
    int isInside(const double3 rx, const double _eps) const;

    string toStr(string spr1="\n------\n", string spr2="\n------\n");
    void print(ostream &output, string spr1="\n------\n", string spr2="\n------\n");
};

class  Surface{
public:
	Triangle *polygons;
	Material *materials;
	int Npol, Nmat;

	Surface(void):
		polygons(NULL), materials(NULL), Npol(0), Nmat(0)
	{}
    Surface(int _n_pol, int _n_mat);

    void resize_clear(int _n_pol, int _n_mat);
    int load_from_file(string surf_filename, string mat_filename);

    void print(ostream &output, string spr1="\n------------------------------\n", string spr2="\n------------------------------\n");
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
string d3ToStr(double3 d, string sp1="(", string sp2=";", string sp3=")");

string error_handl_string, global_logFname = LOG_FILENAME;
