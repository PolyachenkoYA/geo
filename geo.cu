
/**
 * PolyachYA Corporation.  All rights reserved.
 *
 * Please refer to the PolyachYA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "geo.cuh"

long long int n_rays = 0;
long long int n_tot = 0;

// ----------------------------------------- main math-body ------------------------------------------

int rightV(vector<double3>& v, double3 n, double3 v0)
{
	if(v.size() != 4)
		return -1;

	int i,res = -1;
	double cs, max_cs = -2, kv = 1/length(v0);

	for(i = 0; i < 4; ++i){
		if(dot(v0, n)*dot(v[i], n) < 0){ // determine other side of surface
			cs = dot(v[i], v0)/length(v[i])*kv;
			if(max_cs < cs){ // determine reflected ray
				res = i;
				max_cs = cs;
			}
		}
	}

	if(res == -1){
		cerr << "No prover V found (rightV)\n";
		exit(1);
	}
	return res;
}

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy)
{
	vector<double3> v_news;
	v_news.clear();
	v_news.push_back(vx + vy);
	v_news.push_back(vx - vy);
	v_news.push_back(-v_news[1]);
	v_news.push_back(-v_news[0]);
	return v_news[rightV(v_news, n, v_old)];
}

// ----------------------------------------- Other stuff ----------------------------------------------
template<typename T>
void stp(T str){ cout << str << endl; cin.get(); }

template<typename T>
string vectorToStr(vector<T> v, string sp)
{
        string s="";
        for(int i = 0; i < v.size()-1; ++i) s+= (toString(v[i])+sp);
        return s+toString(v[v.size()-1]);
}
template<typename T>
void printVector(ostream &output, vector<T> v, string sp1, string sp2, string sp3)
{
	int sz = v.size();
    output << sp1;
    for(int i = 0; i < sz-1; ++i) output << v[i] << sp2;
    output << v[sz-1] << sp3;
}
template<typename T>
T sumVector(vector<T> v)
{
        T s=0;
        for(int i = 0; i < v.size(); ++i) s+=v[i];
        return s;
}
template <typename T>
string toString(T val)
{
    std::ostringstream oss;
    oss << val;
    return oss.str();
}

template<typename T>
T fromString(const string& s)
{
  std::istringstream iss(s);
  T res;
  iss >> res;
  return res;
}

string toLower(string s)
{
	char d = 'a'-'A';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='A') && (s[i]<='Z')) s[i]+=d;
	return s;
}
string toUpper(string s)
{
	char d = 'A'-'a';
	for(int i = 0; i<s.size(); ++i) if((s[i]>='a') && (s[i]<='z')) s[i]+=d;
	return s;
}
void printD3(ostream &output, double3 r, string sp1, string sp2, string sp3)
{
	output << sp1 << r.x << sp2 << r.y << sp2 << r.z << sp3;
}


// ---------------------------------------- Other math -------------------------------------------------
bool almostEq(double x, double y, double _eps)
{
	return y == 0 ? (abs(x) < _eps) : (abs(x/y - 1) < _eps);
}

vector<double> d3ToV(double3 v){
	double vp[3] = {v.x, v.y, v.z};
	vector<double> vv;
	vv.assign(vp, vp+3);
	return vv;
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
    if(x > _eps){
    	return 1;
    } else if(x < -_eps){
        return -1;
    } else {
        return 0;
    }
}

// ------------------------------------ Triangle ----------------------------------------------
double3 Triangle::getNorm(void) const
{
    return normalize(cross(this->r[1] - this->r[0], this->r[2] - this->r[0]));
}

int Triangle::sg(double3 rx) const
{
    return sgn(dot(rx, this->n));
}

int Triangle::isInside(double3 rx) const
{
    int sg0 = this->sg(cross(this->r[1] - rx, this->r[1] - this->r[0]));
    if(sg0 == 0){
        return 1;
    } else{
        return (sg0 == this->sg(cross(this->r[2] - rx, this->r[2] - this->r[1]))) &&
               (sg0 == this->sg(cross(this->r[0] - rx, this->r[0] - this->r[2])));
    }
}

// -------------------------------------- Ray --------------------------------------------------

void Ray::move(Surf& srf)
{
	n_rays++;
	n_tot++;

	if((this->t + length(this->r)/this->c > Tmax) || (this->A < Amin)){
		n_rays--;
		return;
	}

	//int i;

	// --------------------- pre-geom - find collision point -------------------------
	double3 n, rx;
	pair<double, int> coll_res = this->find_collision(srf);
	double t0 = coll_res.first;
	int i_coll = coll_res.second;
	Triangle trng;

	if(i_coll == -1){
		n_rays--;
		return;
	} else {
		rx = this->r + this->v*t0; // collision point found
		trng = srf.polygons[i_coll]; // collision surface found
		n = trng.n;
	}

	// phi - angle of incidence
	double cos_phi = fabs(dot(n,this->v)/(length(n)*length(this->v)));

	// ------------------------------ phys angles & amplitude ---------------------------
	double sin_phi = cos_sin(cos_phi);
	double p = sin_phi / this->c;
	//double p2 = p*p;
	double sin_p, sin_s, cos_p, cos_s;
	double kp, ks; // As = A*ks, Ap = A*kp

	// snell's law
	// TODO sin>1
	sin_p = p*Cp;
	sin_s = p*Cs;

	// complex math
	// TODO phys
	switch(this->type){
	case PRayType:
		kp = ks = 0.5;
		break;
	case SRayType:
		kp = ks = 0.5;
		break;
	default:
		cerr << "wrong ray type\n";
	}

	// --------------------------------- post-geom - create new rays--------------------------------------
	double3 nr = normalize(cross(n,this->v)); // reflection surface
	double3 v_new;

	// ---- build new P-ray -------
	if(sin_p <= 1){
		cos_p = cos_sin(sin_p);
		v_new = newV(n, this->v, n*cos_p, cross(n, nr)*sin_p);
		Ray p_ray(rx + v_new*eps, v_new, PRayType, this->A*kp, this->t + t0);
		p_ray.move(srf);
	}
	// Ray(const double3 _r, const double3 _ve, const int _type = PRayType, const double _A = 1, const double _t = 0)
	// newV(double3 n, double3 v, double3 vx, double3 vy)

	// ---- build new S-ray -------
	if(sin_s <= 1){
		cos_s = cos_sin(sin_s);
		v_new = newV(n, this->v, n*cos_s, cross(n, nr)*sin_s);
		Ray s_ray(rx + v_new*eps, v_new, SRayType, this->A*ks, this->t + t0);
		s_ray.move(srf);
	}

	n_rays--;
}

pair<double, int> Ray::find_collision(Surf& srf)
{
	int i;
	double t, t_min = Tmax*2;
	Triangle *trngl;
	int i_coll = -1;

	for(i = 0; i < srf.polygons.size(); ++i){ // find collision point

		trngl = &(srf.polygons[i]); // so we don't have to call [i] every time, also it's shorter
		t = dot(trngl->r[0] - this->r, trngl->n) / dot(this->v, trngl->n); // find time of collision

		if((0 < t) && (t < t_min)){ // if the possible collision can happen (t>0) and if better than the one we already have(t<t_min)
			if(trngl->isInside(this->r + this->v*t)){ // if it's really the point, then save it
				t_min = t;
				i_coll = i;
			}
		}
	}

	return make_pair(t_min, i_coll);
}

void Ray::print(ostream &output, string spr)
{
	output << spr;
	output << (this->type == PRayType ? "P-type" : "S-type") << "\n";
	output << "c = " << this->c << "; A = " << this->A << "; t = " << this->t << "\n";
	output << "r = ";
	printD3(output, this->r);
	output << "v = ";
	printD3(output, this->v);
	output << spr;
}

