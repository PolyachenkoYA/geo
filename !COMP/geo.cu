
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
		if(dot(v0, n)*dot(v[i], n) < 0){ // exclude other side of surface
			cs = dot(v[i], v0)/length(v[i])*kv;
			if(cs > max_cs){ // determine reflected ray
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
char bool2int(bool b){ return b ? 1 : 0; }

template<typename T>
//int stp(T str)
void stp(T str)
{
	cerr << str << endl;
	cin.get();
	//int r;
	//cin >> r;
	//return r;
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
	output << d3ToStr(r, sp1, sp2, sp3);
}

template<typename T>
int CHECK(int n, T s)
{
    if(!n){ return 0; }
    ofstream Fout(global_logFname.c_str(), ios::app);
    if(!Fout){ return CANT_OPEN_FILE_FOR_WRITING; }

    if(n != SAY_IT){
    	string _s = "error #" + toString(n) + "\n   message:\n" + error_handl_string + "\n";
    	Fout << "\n" << _s;
    	cout << _s;
    }
    switch(n){
    	case SAY_IT: Fout << s; break;

    	case CANT_OPEN_FILE_FOR_READING: Fout << "Can't open file\n" << s << "\nfor reading"; break;
    	case CANT_OPEN_FILE_FOR_WRITING: Fout << "Can't open file\n" << s << "\nfor writing"; break;

    	case WRONG_SET_OF_POSSIBLE_REFLECTED_RAYS: Fout << "Wrong set of possible reflected rays (v.size != 4)\n"; break;
    	case NO_REFLECTED_RAY_FOUND: Fout << "No reflected ray found\n"; break;
    	case WRONG_RAY_TYPE: Fout << "Wrong ray type: " << s; break;

        default: Fout << "Unknown error #" << n << "\nerror message:\n" << s;
    }
    if(n != SAY_IT){ Fout << "\n"; }
    Fout.close();

    return n;
}

template<typename T>
void SAY_LOG(T s)
{
	CHECK(SAY_IT, s);
}

// ---------------------------------------- Other math -------------------------------------------------
bool almostEq(double x, double y, double _eps)
{
	return y == 0 ? (abs(x) < _eps) : (abs(x/y - 1) < _eps);
}

string d3ToStr(double3 d, string sp1, string sp2, string sp3)
{
	return sp1 + toString(d.x) + sp2 + toString(d.y) + sp2 + toString(d.z) + sp3;
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
    if(x >= _eps){
    	return 1;
    } else if(x <= -_eps){
        return -1;
    } else {
        return 0;
    }
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Material ----------------------------------------------
// --------------------------------------------------------------------------------------------

bool Material::isEq(const Material* m2, const double _eps) const
{
	return m2 == NULL ? 0 : (almostEq(this->Cp, m2->Cp, _eps) && almostEq(this->Cs, m2->Cs, _eps));
}

string Material::toStr(string spr1, string spr2)
{
	return spr1 + "Cp = " + toString(this->Cp) + "; Cs = " + toString(this->Cs) + spr2;
}

void Material::print(ostream &output, string spr1, string spr2)
{
	output << this->toStr(spr1, spr2);
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Triangle ----------------------------------------------
// --------------------------------------------------------------------------------------------
Triangle::Triangle(){
	this->clear_params();
}
Triangle::Triangle(const double3 r1, const double3 r2, const double3 r3){
	this->clear_params();
    this->r[0] = r1;
    this->r[1] = r2;
    this->r[2] = r3;
    this->n = this->getNorm();
}
Triangle::Triangle(const double3* r_new)
{
	this->clear_params();
	int i;
	for(i = 0; i < 3; ++i){
		this->r[i] = r_new[i];
	}
	this->n = this->getNorm();
}
Triangle::Triangle(const double3 *vertex, const int* ind)
{
	this->clear_params();
	int i;
	for(i = 0; i < 3; ++i){
		this->i_vrtx[i] = ind[i];
	}
	for(i = 0; i < 3; ++i){
		this->r[i] = vertex[this->i_vrtx[i]];
	}
	this->n = this->getNorm();
}
Triangle::Triangle(const double3 *vertex, const int i0, const int i1, const int i2)
{
	this->clear_params();
	this->i_vrtx[0] = i0;
	this->i_vrtx[1] = i1;
	this->i_vrtx[2] = i2;

	for(int i = 0; i < 3; ++i){
		this->r[i] = vertex[this->i_vrtx[i]];
	}
	this->n = this->getNorm();
}

void Triangle::clear_params()
{
	this->mat[0] = this->mat[1] = NULL;
	this->i_vrtx[0] = this->i_vrtx[1] = this->i_vrtx[2] = 0;
	this->r[0] = this->r[1] = this->r[2] = this->n = make_double3(0,0,0);
}

double3 Triangle::getNorm(void) const
{
    return normalize(cross(this->r[1] - this->r[0], this->r[2] - this->r[0]));
}

int Triangle::sg(const double3 rx, const double _eps) const
{
    return sgn(dot(rx, this->n), _eps);
}

int Triangle::isInside(const double3 rx, const double _eps) const
{
    int sg0 = this->sg(cross(this->r[1] - rx, this->r[1] - this->r[0]), _eps);
    return sg0 == 0 ? 0 : ((sg0 == this->sg(cross(this->r[2] - rx, this->r[2] - this->r[1]), _eps)) &&
    		              (sg0 == this->sg(cross(this->r[0] - rx, this->r[0] - this->r[2]), _eps)));
}

string Triangle::toStr(string spr1, string spr2)
{
	return spr1 + "vertices coords are:\n" +
		   d3ToStr(this->r[0]) + "\n" +
		   d3ToStr(this->r[1]) + "\n" +
		   d3ToStr(this->r[2]) + "\n" +
		   "i_vrtx = " + toString(this->i_vrtx[0]+1) + ";" + toString(this->i_vrtx[1]+1) + ";" + toString(this->i_vrtx[2]+1) + "\n" +
		   "n = " + d3ToStr(this->n) + "\n" +
		   "mat0 : " + this->mat[0]->toStr() + "\n"
		   "mat1 : " + this->mat[1]->toStr() + spr2;
}

void Triangle::print(ostream &output, string spr1, string spr2)
{
	output << this->toStr(spr1, spr2);
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Ray ---------------------------------------------------
// --------------------------------------------------------------------------------------------

void Ray::move(Surface* srf, Params *prm)
{
	this->print(cout);
	prm->n_alive_rays++;
	prm->n_total_rays++;

	if((this->t + length(this->r)/this->c > prm->Tmax) || (this->A < prm->Amin)){
		prm->n_alive_rays--;
		return;
	}

	//int i;

	// --------------------- pre-geom - find collision point -------------------------
	double3 n, rx;
	pair<double, int> coll_res = this->find_collision(srf, prm);
	// pair<double, int> Ray::find_collision(Surface* srf, Params* prm)
	double t0 = coll_res.first;
	int i_coll = coll_res.second;
	Triangle *trng;

	if(i_coll == -1){
		prm->n_alive_rays--;
		return;
	} else {
		rx = this->r + this->v*t0; // collision point found
		trng = &(srf->polygons[i_coll]); // collision surface found
		n = trng->n;
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
	Material *mat_from = trng->mat[ bool2int(trng->sg(trng->r[0] - this->r, 0) < 0) ];
	sin_p = p * mat_from->Cp;
	sin_s = p * mat_from->Cs;

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
		error_handl_string = "Wrong ray type: " + toString(this->type) + this->toStr();
		CHECK(WRONG_RAY_TYPE, this->type);
		exit(1);
	}

	// --------------------------------- post-geom - create new rays--------------------------------------
	double3 nr = normalize(cross(n,this->v)); // reflection surface
	double3 v_new;

	// ---- build new P-ray -------
	if(sin_p <= 1){
		cos_p = cos_sin(sin_p);
		v_new = newV(n, this->v, n*cos_p, cross(n, nr)*sin_p) * mat_from->Cp;
		Ray p_ray(rx + v_new * prm->eps, v_new, PRayType, this->A*kp, this->t + t0);
		p_ray.move(srf, prm);
	}
	// Ray(const double3 _r, const double3 _v, const int _type = BaseRayType, const double _A = 1, const double _t = 0):
	// newV(double3 n, double3 v, double3 vx, double3 vy)

	// ---- build new S-ray -------
	if(sin_s <= 1){
		cos_s = cos_sin(sin_s);
		v_new = newV(n, this->v, n*cos_s, cross(n, nr)*sin_s) * mat_from->Cs;
		Ray s_ray(rx + v_new*prm->eps, v_new, SRayType, this->A*ks, this->t + t0);
		s_ray.move(srf, prm);
	}

	prm->n_alive_rays--;
}

pair<double, int> Ray::find_collision(Surface* srf, Params* prm)
{
	int i;
	double t, t_min = prm->Tmax * 2;
	Triangle *trngl;
	int i_coll = -1;

	for(i = 0; i < srf->Npol; ++i){ // find collision point

		trngl = &(srf->polygons[i]); // so we don't have to call [i] every time, also it's shorter
		t = dot(trngl->r[0] - this->r, trngl->n) / dot(this->v, trngl->n); // find time of collision

		if((SYS_EPS < t) && (t < t_min)){
		// if the possible collision can happen (t>0) and if it's better than the one we already have (t<t_min)
			if(trngl->isInside(this->r + this->v*t, prm->eps)){ // if it's really the point, then save it
				t_min = t;
				i_coll = i;
			}
		}
	}

	return make_pair(t_min, i_coll);
}

string Ray::toStr(const string spr)
{
	return spr + (this->type == PRayType ? "P-type" : "S-type") +
			"\nc = " + toString(this->c) + "; A = " + toString(this->A) + "; t = " + toString(this->t) +
			"\nr = " + d3ToStr(this->r) +
			"\nv = " + d3ToStr(this->v) + spr;
}

void Ray::print(ostream &output, string spr)
{
	output << this->toStr(spr);
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Surface -----------------------------------------------
// --------------------------------------------------------------------------------------------

Surface::Surface(int _n_pol, int _n_mat)
{
	this->resize_clear(_n_pol, _n_mat);
}

void Surface::resize_clear(int _n_pol, int _n_mat)
{
	if(_n_pol != this->Npol){
		delete[] this->polygons;
		this->polygons = (_n_pol == 0 ? NULL : (new Triangle[_n_pol]));
		this->Npol = _n_pol;
	}
	if(_n_mat != this->Nmat){
		delete[] this->materials;
		this->materials = (_n_mat == 0 ? NULL : (new Material[_n_mat]));
		this->Nmat = _n_mat;
	}
}

void Surface::print(ostream &output, string spr1, string spr2)
{
	int i;
	output << spr1 << "Nmat = " << this->Nmat << "\n";
	for(i = 0; i < this->Nmat; ++i){
		this->materials[i].print(output, toString(i+1) + ") ");
	}
	output << "Npol = " << this->Npol << "\n";
	for(i = 0; i < this->Npol; ++i){
		this->polygons[i].print(output, "\n------------\n" + toString(i+1) + "\n");
	}
	output << spr2;
}

int Surface::load_from_file(string surf_filename, string mat_filename)
{
	ifstream input;
	int i, i2, i3, i4, i5;
	string buf_str;

	// --------------------------------------- read materials -----------------------------------------------
	input.open(mat_filename);
	if(!input){
		input.close();
		error_handl_string = "material-file " + mat_filename + " is missing\n";
		return CANT_OPEN_FILE_FOR_READING;
	}

	int _n_mat;
	input >> _n_mat;
	if(_n_mat < 2){
		error_handl_string = "Nmat = " + toString(_n_mat) + "; must be at least 2\n";
		return LESS_2_MATERIALS;
	}
	this->resize_clear(this->Npol, _n_mat);
	for(i = 0; i < this->Nmat; ++i){
		input >> this->materials[i].Cs >> this->materials[i].Cp;
	}
	input.close();


	// ---------------------------------------------------------------------------------------------------
	// ------------------------------------- read surface ------------------------------------------------
	// ---------------------------------------------------------------------------------------------------
	input.open(surf_filename);
	if(!input){
		input.close();
		error_handl_string = "surface-file " + surf_filename + " is missing\n";
		return CANT_OPEN_FILE_FOR_READING;
	}
	double3 *vertex;
	int Nvertex;
	// --------------------------------------- read vertex -----------------------------------------------
	input >> Nvertex;
	vertex = new double3[Nvertex];
	for(i = 0; i < Nvertex; ++i){
		input >> vertex[i].x >> vertex[i].y >> vertex[i].z;
	}

	// ------------------------------------- read tetrahedrons & triangles ---------------------------------------------
	int Ntet, i_mat, i_vertex[4];

	input >> Ntet;

	vector<Triangle> read_polygons;
	int n_overlap_vrtx, i_curr_mat;
	int tet_surf_ind[4][3] =
	// tet_surf_ind[i][:] - all vertices except i-th one -
	// the indices of the i-th triangle-sub-surf of the abstract tetrahedron
	{
	  {1, 2, 3},
	  {0, 2, 3},
	  {0, 1, 3},
	  {0, 1, 2}
	};

	for(i = 0; i < Ntet; ++i){
		for(i2 = 0; i2 < 4; ++i2){ // read vertices of tetrahedron_i
			input >> i_vertex[i2];
			--i_vertex[i2]; // arrays from 0, in file-format from 1
		}
		input >> i_mat; // read material of tetrahedron_i

		for(i2 = 0; i2 < 4; ++i2){           // for each tetrahedron sub-surface (tet_surf_ind[i2][:]) check overlaps
			for(i3 = 0; i3 < read_polygons.size(); ++i3){       // with all existing ones
				n_overlap_vrtx = 0;
				for(i4 = 0; i4 < 3; ++i4){
					for(i5 = 0; i5 < 3; ++i5){
						/*
						 * i - index of the current tetrahedron
						 * i2 - index of surface of current tetrahedron
						 * i3 - for all existing sufraces - index of some existing surface to check a pair tetr_surf<->existing_surf
						 * i4 - index of the vertex of the current surf of the current tetr
						 * i5 - index of the vertex of the current existing triangle
						 */
						if(i_vertex[ tet_surf_ind[i2][i4] ] == read_polygons[i3].i_vrtx[i5]){
							++n_overlap_vrtx;
							break;
						}
					}
				}
				if(n_overlap_vrtx >= 3){
					break;
					/*
					 * current (i2) triangle-sub-surf can't overlap with more than 1 existing triangle.
					 * It's found, so no overlapping with all further ones.
					 * */
				}
			}
			if(n_overlap_vrtx > 3){
				error_handl_string = "n_overlap_vrtx = " + toString(n_overlap_vrtx) + "\n i = " + toString(i) + "; i2 = " + toString(i2) + "; i3 = " + toString(i3);
				delete[] vertex;
				read_polygons.clear();
				input.close();
				return TOO_MANY_OVERLAPED_VERTICES;
			}
			if(i3 == read_polygons.size()){ // current (i2) triangle-sub-surf is new to the all existing ones
				read_polygons.push_back(Triangle(vertex, i_vertex[ tet_surf_ind[i2][0] ],
						                                 i_vertex[ tet_surf_ind[i2][1] ],
						                                 i_vertex[ tet_surf_ind[i2][2] ]));
			}
			// else means already used polygons similar to the current (i2) one was found
			// read_polygons[i3].i_vrtx[:] == i_vertex[ tet_surf_ind[i2][:] ]

			/*
			 * if the current triangle is a new one, the i3 is the number of polygons BEFORE new one was added,
			 * so after all i3 is the index of the last added triangle. In this case i3 == read_polygons.size()-1
			 *
			 * otherwise (existing similar one was found) i3 is the index of the triangle identical to the current one.
			 * In this case we add 2nd matetrial to the existing polygon and don't add a new triangle to polygons
			 *
			 * So in both cases i3 is the index of the triangle to add the current material to
			 */
			i_curr_mat = bool2int(read_polygons[i3].sg(vertex[ i_vertex[i2] ] - read_polygons[i3].r[0], 0) > 0);
			// i2 is the missing index in tet_surf_ind[i2][:] set, so i_vertex[i2] is the 4th vertex of the current tetrahedron
			read_polygons[i3].mat[i_curr_mat] = this->materials + i_mat;

			if(read_polygons[i3].mat[i_curr_mat]->isEq(read_polygons[i3].mat[bool2int(!i_curr_mat)])){
				// read_polygons[i3] is a formal border with no real meaning
				read_polygons.erase(read_polygons.begin()+i3);
			}
		}
	}
	/*
	 * As a result we have read_polygons with all true-surface triangles of the system
	 * Now we need to add background material and copy read polygons to the main array
	 */
	for(i = 0; i < read_polygons.size(); ++i){ // if some materials are unset, then set them to background material
		for(i2 = 0; i2 < 2; ++i2){
			if(read_polygons[i].mat[i2] == NULL){
				read_polygons[i].mat[i2] = this->materials; // &(material[0])
			}
		}
		if(read_polygons[i].mat[0]->isEq(read_polygons[i].mat[1])){
			// read_polygons[i] is a formal border with no real meaning
			// actually it can't happen here if everything before worked fine, but lets leave it this way for now
			// TODO error polygon with no material
			read_polygons.erase(read_polygons.begin()+i);
		}
	}

	this->resize_clear(read_polygons.size(), this->Nmat);
	for(i = 0; i < this->Npol; ++i){ // copy read materials
		this->polygons[i] = read_polygons[i];
	}

	delete[] vertex;
	read_polygons.clear();
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

	string buf_s;
	std::getline(input, buf_s); // read comment line

	double _eps;
	input >> this->Nrays >> this->Tmax >> this->Amin >> _eps;
	this->eps = _eps > 0 ? _eps : SYS_EPS;
	input.close();

	return 0;
}

void Params::print(ostream &output, string spr)
{
	for(int i = 0; i < this->paramFHead.size(); ++i){
		output << this->paramFHead[i] << " | ";
	}
	output << spr << this->Nrays << spr
		   << this->Tmax << spr
		   << this->Amin << spr
		   << this->eps << spr
		   << "model_name: |" << this->model_name << "|" << spr
		   << "alive_rays: " << this->n_alive_rays << "; total_rays: " << this->n_total_rays << spr;
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
