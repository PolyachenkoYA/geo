
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

string error_handl_string;

// ----------------------------------------- main math ------------------------------------------

int rightV(vector<double3>& v, double3 n, double3 v0)
{
	if(v.size() != 4){
		CHECK(WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS, 0);
		exit(1);
	}


	int i,res = -1;
	double cs, max_cs = -2, kv = 1/length(v0);

	// reflected ray have the least angle with v0 |=> the beggest cos(v_new,v0)
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
		CHECK(NO_REFLECTED_RAY_FOUND, 0);
		exit(1);
	}
	return res;
}

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy)
{
	vector<double3> v_news;
	v_news.clear();
	// we don't know the direction of n and nr, so we have to check all 4 possible variants
	// to choose the one that really is the physicaly reflected ray
	v_news.push_back(vx + vy);
	v_news.push_back(vx - vy);
	v_news.push_back(-v_news[1]);
	v_news.push_back(-v_news[0]);
	return v_news[rightV(v_news, n, v_old)];
}

// ----------------------------------------- Other stuff ----------------------------------------------
template<typename T>
void stp(T str)
{
	cout << str << endl; cin.get();
}

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

int find_partition(ifstream &input, string mask, bool whole_file = 0)
{
	if(whole_file){
		input.clear();
		input.seekg(0, ios::beg);
	}
	int res = 0;

	char buf_chars[256];
	string buf_str;
	do{
		++res;
		input.getline(buf_chars, 255);
		buf_str = string(buf_chars);
		// getline doesn't remove '\n' in the end of a line, so we compare mask with buf_str[0:-1]
	}while(buf_str.compare(0, buf_str.size()-1, mask) && !input.eof());

	if(buf_str.compare(0, buf_str.size()-1, mask)){
		CHECK(WRONG_MSH_FILE_ID, "No |" + mask + "| partition found in msh file");
		exit(1);
	}

	return res;
}

template<typename T>
int CHECK(int n, T s, string logFname)
{
    if(!n){ return 0; }
    ofstream Fout(logFname.c_str(), ios::app);
    if(!Fout){ return CANT_OPEN_FILE_FOR_WRITING; }

    if(n != SAY_IT){
    	Fout << "\nerror #" << toString(n) << "\n";
    	cout << "error #" << toString(n) << "\n";
    }
    switch(n){
    	case SAY_IT: Fout << s; break;

    	case CANT_OPEN_FILE_FOR_READING: Fout << "Can't open file\n" << s << "\nfor reading"; break;
    	case CANT_OPEN_FILE_FOR_WRITING: Fout << "Can't open file\n" << s << "\nfor writing"; break;

    	case TOO_BIG_NODE_IND: Fout << "Too big node ID in file\n" << s << "\n" << error_handl_string; break;
    	case TOO_BIG_ELEMENT_IND: Fout << "Too big element ID in file\n" << s << "\n" << error_handl_string; break;
    	case WRONG_MSH_FILE_ID: Fout << s; break;

    	case WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS: Fout << "Wrong set of possible reflected rays (v.size != 4)\n"; break;
    	case NO_REFLECTED_RAY_FOUND: Fout << "No reflected ray found\n"; break;

        default: Fout << "Unknown error #" << n;
    }
    if(n != SAY_IT){ Fout << "\n"; }
    Fout.close();

    return n;
}

template<typename T>
void SAY_LOG(T s, string logFname)
{
	CHECK(SAY_IT, s, logFname);
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

// --------------------------------------------------------------------------------------------
// ------------------------------------ Triangle ----------------------------------------------
// --------------------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------------------
// ------------------------------------ Ray ---------------------------------------------------
// --------------------------------------------------------------------------------------------

void Ray::move(Surface& srf)
{
	this->print(cout);
	prm.n_alive_rays++;
	prm.n_total_rays++;

	if((this->t + length(this->r)/this->c > prm.Tmax) || (this->A < prm.Amin)){
		prm.n_alive_rays--;
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
		prm.n_alive_rays--;
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
	sin_p = p * prm.Cp;
	sin_s = p * prm.Cs;

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
		Ray p_ray(rx + v_new*prm.eps, v_new, PRayType, this->A*kp, this->t + t0);
		p_ray.move(srf);
	}
	// Ray(const double3 _r, const double3 _ve, const int _type = PRayType, const double _A = 1, const double _t = 0)
	// newV(double3 n, double3 v, double3 vx, double3 vy)

	// ---- build new S-ray -------
	if(sin_s <= 1){
		cos_s = cos_sin(sin_s);
		v_new = newV(n, this->v, n*cos_s, cross(n, nr)*sin_s);
		Ray s_ray(rx + v_new*prm.eps, v_new, SRayType, this->A*ks, this->t + t0);
		s_ray.move(srf);
	}

	prm.n_alive_rays--;
}

pair<double, int> Ray::find_collision(Surface& srf)
{
	int i;
	double t, t_min = prm.Tmax*2;
	Triangle *trngl;
	int i_coll = -1;

	for(i = 0; i < srf.size; ++i){ // find collision point

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

// --------------------------------------------------------------------------------------------
// ------------------------------------ Surface -----------------------------------------------
// --------------------------------------------------------------------------------------------

Surface::Surface(int _n)
{
	if(_n == 0){
		this->polygons = NULL;
		this->size = 0;
	} else {
		this->resize_clear(_n);
	}
}

void Surface::resize_clear(int _n)
{
	if(this->polygons != NULL){
		delete[] this->polygons;
	}
	if(_n == 0)
		return;

	this->polygons = new Triangle[_n];
	this->size = _n;
}

int Surface::load_from_file(string filename)
{
	ifstream input(filename);
	if(!input){
		input.close();
		return CANT_OPEN_FILE_FOR_READING;
	}

	int i, err_handl;
	double buf_d;
	string buf_str;
	double3 *nodes;
	int Nnodes, k;

	int node_partition_i = find_partition(input, "$Nodes", 0);
	input >> Nnodes;
	nodes = new double3[Nnodes];
	for(i = 0; i < Nnodes; ++i){ // load nodes
		input >> k;
		--k; // k begins from 1, array begins from 0
		// actually the .gmh format says that nodes don't have to be ordered, so it's possible for 'k' to be > Nnodes.
		// But for our work lets state (at least for now) that it can't be so.
		// All nodes [1;N] have to be set before N+1 node will be set
		if(k >= Nnodes){
			delete[] nodes;
			input.close();
			error_handl_string = "local_line_ind = " + toString(i) + "; k = " + toString(k+1);
			return TOO_BIG_NODE_IND;
		}
		input >> nodes[k].x >> nodes[k].y >> nodes[k].z;
	}

	int elements_partition_i = find_partition(input, "$Elements", 0);
	int Nel, elType, k_r1, k_r2, k_r3;
	input >> Nel;
	this->resize_clear(Nel);
	for(i = 0; i < this->size; ++i){ // build elements using already read nodes
		input >> k;
		--k; // k begins from 1, array begins from 0
		if(k >= this->size){
			delete[] nodes;
			input.close();
			error_handl_string = "local_line_ind = " + toString(i) + "; k = " + toString(k+1);
			return TOO_BIG_ELEMENT_IND;
		}

		input >> elType;

		if(elType == 2){
			input >> buf_d; // teg_number
			input >> buf_d; // 1st teg - number of the physical entity to which the element belongs
			input >> buf_d; // 2nd teg - number of the elementary geometrical entity to which the element belongs

			input >> k_r1 >> k_r2 >> k_r3;
			// k begins from 1, array begins from 0
			--k_r1;
			--k_r2;
			--k_r3;
			if((k_r1 >= Nnodes) || (k_r2 >= Nnodes) || (k_r3 >= Nnodes)){
				delete[] nodes;
				input.close();
				error_handl_string = "local_line_ind = " + toString(i) + "; k1 = " + toString(k_r1+1) + "; k2 = " + toString(k_r2+1) + "; k3 = " + toString(k_r3+1);
				return TOO_BIG_NODE_IND;
			}

			this->polygons[k] = Triangle(nodes[k_r1], nodes[k_r2], nodes[k_r3]);
		} else {
			// kind of error, but ok, let's just warn the user about not complitely supported format
			buf_str = "line " + toString(i) + " in file\n" + filename + "\nis't a triangle\nOnly triangles are supported yet.\n";
			SAY_LOG(buf_str);
			cout << buf_str;
		}
	}

	delete[] nodes;
	input.close();
	return 0;
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Surface -----------------------------------------------
// --------------------------------------------------------------------------------------------


int Params::load_from_file(string filename)
{
	ifstream input(filename);
	if(!input){
		input.close();
		return CANT_OPEN_FILE_FOR_READING;
	}

	double _eps;
	input >> this->Nrays >> this->Tmax >> this->Amin >> this->Cp >> this->Cs >> _eps;
	this->eps = _eps > 0 ? _eps : SYS_EPS;

	if(!input.eof()){
		char buf_str[256];
		input.getline(buf_str, 255);
		this->msh_filename = string(buf_str);
	}

	input.close();
	return 0;
}

int Params::print(ostream &output, string spr)
{
	output << this->Nrays << spr
		   << this->Tmax << spr
		   << this->Amin << spr
		   << this->Cp << spr
		   << this->Cs << spr
		   << this->eps << spr
		   << this->msh_filename << spr;

	return 0;
}

int Params::save_to_file(string filename)
{
	ofstream output(filename);
	if(!output){
		return CANT_OPEN_FILE_FOR_WRITING;
	}

	this->print(output);

	output.close();
	return 0;
}
