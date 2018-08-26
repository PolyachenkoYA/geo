#pragma once

/*
 * Space.cuh
 *
 *  Created on: Jul 7, 2016
 *      Author: polyachyadt
 */

//#ifndef GEO_CUH_
//#define GEO_CUH_

#include "const.cuh"
#include "general_math.cuh"
#include "spec_math.cuh"
#include "format.cuh"

class Ray;
class Triangle;
class Surface;
class Params;
class Material;
class RegisteredRay;
class Detector;

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
	double Tmax, Amin, dt, f, tau, B;
	int Nrays, Nfrm;
	long long int n_alive_rays, n_total_rays;
	bool use_det, draw_mov;
	int prnt_mode;

	double3 Xmax, Xmin, dX;
	int3 Nslc;

	string model_name, tet_filename, prm_filename, material_filename, res_filename;

	vector<string> paramFHead;

	Params(void):
		eps(SYS_EPS), Tmax(0), dt(0), Amin(0), f(0), tau(0), Nrays(0), n_alive_rays(0), n_total_rays(0), model_name("sys_tst")
	{
		this->paramFHead = {"Nrays","Tmax","dt","Amin","f","tau","eps","use_det","draw_mov","prnt_mode",
				            "Xmin.x","Xmin.y","Xmin.z","Xmax.x","Xmax.y","Xmax.z","dX.x","dX.y","dX.z"};
	}
	~Params(void){}

	int load_from_file(string filename);
	void print(ostream &output);
	void print_full(ostream &output, string spr);
	int save_to_file(string filename);
};

class Ray{
public:
    int type;
    double c, A, t;
    double3 r, v, polar;

    Ray(void):
            type(PRayType), c(0), A(0), t(0)
    { v = r = polar = make_double3(0,0,0); }
    Ray(Ray *_r);
    Ray(const double3 _r, const double3 _v, const double3 _polar, const int _type = PRayType, const double _A = 1, const double _t = 0);
    ~Ray(){}

    void print(ostream &output, string spr = "\n-----\n");
    string ToStr(const string spr = "\n-----\n");
    int move(Surface* srf, Params *prm);
    pair<double, int> find_collision(Surface* srf, Params* prm);
    RegisteredRay toRegRay(void);
    void add(Ray *ray2);
};

class Triangle{
public:
    double3 r[3];
    int i_vrtx[3]; // r[i] == nodes[i_vrtx[i]]
    double3 n; // n = [(r1-r0)x(r2-r0)]

    Material* mat[2];
    // mat[1] is in the side where normal-vector points, mat[0] is in the other (back) side

    bool is_absorber;
    Detector *detector;
    // don't free these pointers because they are just pointers, not a unique arrays or smth like that

    Triangle();
    Triangle(const double3 r1, const double3 r2, const double3 r3);
    Triangle(const double3* r_new);
    Triangle(const double3 *vertex, const int* ind);
    Triangle(const double3 *vertex, const int i0, const int i1, const int i2);

    void clear_params();
    double3 getNorm(void) const;
    int sg(const double3 rx, const double _eps = SYS_EPS) const;
    int isInside(const double3 rx, const double _eps) const;

    string ToStr(string spr1="\n------\n", string spr2="\n------\n");
    void print(ostream &output, string spr1="\n------\n", string spr2="\n------\n");
};

class RegisteredRay{
public:
	int type;
	double c, A, t;

	RegisteredRay(int _type = 0, double _A = 0, double _t = 0, double _c = 0):
		type(_type), A(_A), t(_t), c(_c)
	{}
	RegisteredRay(Ray *_ray){
		this->type = _ray->type;
		this->A = _ray->A;
		this->t = _ray->t;
		this->c = _ray->c;
	}
};

class Frame{
public:
	vector<Ray> regRays;

	Frame(void){
		this->regRays.clear();
	}

	int saveToFile(ostream &output, Params *prm);
};

class Detector{
public:
	vector<RegisteredRay> regRays;

	Detector(void){
		this->regRays.clear();
	}

	void saveInfo(ostream &output, Params *prm);
	double regValue(double t0, Params *prm);
	double peakFnc(double x, double B, double A = 1);
};

class  Surface{
public:
	Triangle *polygons;
	Material *materials;
	Detector *detectors;
	Frame *frames; // TODO Nfrm here
	int Npol, Nmat, Ndet;
	double absorbedP, absorbedS, totalE0;

	Surface(void):
		polygons(NULL), materials(NULL), detectors(NULL), frames(NULL), Npol(0), Nmat(0), Ndet(0)
	{}
    Surface(int _n_pol, int _n_mat, int _n_det);
    ~Surface(void){
    	delete[] this->polygons;
    	delete[] this->materials;
    	delete[] this->detectors;
    	delete[] this->frames;
    }

    void resize_clear(int _n_pol, int _n_mat, int _n_det);
    int load_from_file(string surf_filename, string mat_filename, Params *prm);
    Triangle* findPolygon(int* i_vertex);

    int saveDetectorInfo(string filename, Params *prm);
    int saveMovie(Params *prm);
    void print(ostream &output, string spr1="\n------------------------------\n", string spr2="\n------------------------------\n");
};
