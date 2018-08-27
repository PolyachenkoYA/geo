
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

// --------------------------------------------------------------------------------------------
// ------------------------------------ Material ----------------------------------------------
// --------------------------------------------------------------------------------------------

bool Material::isEq(const Material* m2, const double _eps) const
{
	return m2 ? (almostEq(this->Cp, m2->Cp, _eps) && almostEq(this->Cs, m2->Cs, _eps)) : 0;
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
	this->mat[0] = this->mat[1] = nullptr;
	this->detector = nullptr;
	this->is_absorber = 0;
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

string Triangle::ToStr(string spr1, string spr2)
{
	return spr1 + "vertices coords are:\n" +
		   toStr(this->r[0]) + "\n" +
		   toStr(this->r[1]) + "\n" +
		   toStr(this->r[2]) + "\n" +
		   "i_vrtx = " + toString(this->i_vrtx[0]+1) + ";" + toString(this->i_vrtx[1]+1) + ";" + toString(this->i_vrtx[2]+1) + "\n" +
		   "n = " + toStr(this->n) + "\n" +
		   "mat0 : " + this->mat[0]->toStr() + "\n" +
		   "mat1 : " + this->mat[1]->toStr() +  "\n" +
		   "absorber: " + (this->is_absorber ? "yes" : "no") + "\n" +
		   "detector: " + (this->detector ? "yes" : "no") + spr2;
}

void Triangle::print(ostream &output, string spr1, string spr2)
{
	output << this->ToStr(spr1, spr2);
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Ray ---------------------------------------------------
// --------------------------------------------------------------------------------------------

Ray::Ray(Ray *_r)
{
	this->type = _r->type;
	this->c = _r->c;
	this->A = _r->A;
	this->t = _r->t;
	this->r = _r->r;
	this->v = _r->v;
	this->polar = _r->polar;
	this->next = _r->next;
}
Ray::Ray(const double3 _r, const double3 _v, const double3 _polar, const int _type, const double _A, const double _t, Ray *_next):
        type(_type), A(_A), t(_t), next(_next)
{
	this->r = _r;
    this->v = _v;
    this->polar = _polar;
    this->c = length(this->v);
}

RegisteredRay Ray::toRegRay(void)
{
	return RegisteredRay(this);
}

void Ray::add(Ray *ray2)
{
	this->v = (this->v*this->A + ray2->v*ray2->A) / (this->A + ray2->A);
	this->A += ray2->A;
	this->c = ray2->c;
	this->type = ray2->type;
	this->polar = ray2->polar;
	// t & r were determined before
}

int Ray::move(Surface* srf, Params *prm, RaysFront *rays)
{
	if((this->t > prm->Tmax*(1 + prm->eps)) || (this->A < prm->Amin)){
		rays->quit_ray();
		return 0;
	}

	int i;

	// --------------------- pre-geom - find collision point -------------------------
	double3 n, rx;
	pair<double, int> coll_res = this->find_collision(srf, prm);
	// pair<double, int> Ray::find_collision(Surface* srf, Params* prm)
	double dlt_t = coll_res.first;
	int i_coll = coll_res.second;
	Triangle *trng;

	if(i_coll == -1){ // no collision found, so the ray just runs away from the surface
		rays->quit_ray();
		return 0;
	} else {
		trng = &(srf->polygons[i_coll]); // collision surface found
		double3 dlt_r = this->v*dlt_t;

		if(prm->draw_mov){ // draw frames
			double t_new = this->t + dlt_t;
			i = int(this->t / prm->dt) + 1;
			double frame_t = i * prm->dt;
			int3 Xi;
			int ind;
			double3 r_curr;

			// recreate positions of the ray between collisions
			while((frame_t < t_new) && (i < prm->Nfrm)){
				r_curr = this->r + dlt_r * ((frame_t - this->t) / dlt_t); // r when t = frame_t
				// prevent rays from runnig out of the whole system in case of
				// __________|_______|____________________|_____
				//          t_i    t_new                 t_i+1
				if(r_curr.x >= prm->Xmax.x) r_curr.x = prm->Xmax.x - prm->eps;
				if(r_curr.y >= prm->Xmax.y) r_curr.y = prm->Xmax.y - prm->eps;
				if(r_curr.z >= prm->Xmax.z) r_curr.z = prm->Xmax.z - prm->eps;
				if(r_curr.x <= prm->Xmin.x) r_curr.x = prm->Xmin.x + prm->eps;
				if(r_curr.y <= prm->Xmin.y) r_curr.y = prm->Xmin.y + prm->eps;
				if(r_curr.z <= prm->Xmin.z) r_curr.z = prm->Xmin.z + prm->eps;

				switch(prm->prnt_mode){
				case RAW_DATA_MODE:
					srf->frames[i].regRays.push_back(Ray(r_curr, this->v, this->polar, this->type, this->A, frame_t));
					// Ray(const double3 _r, const double3 _v, const int _type = BaseRayType, const double _A = 1, const double _t = 0)
					break;
				case FRAMES_DATA_MODE:
					Xi = get_Nslc(r_curr - prm->Xmin, prm->dX);
					if((Xi.x >= prm->Nslc.x) || (Xi.y >= prm->Nslc.y) || (Xi.z >= prm->Nslc.z)){
						error_handl_string = "r_curr = " + toStr(r_curr) + "; Xi = " + toStr(Xi) + "; Nslc = " + toStr(prm->Nslc) + "; t = " + toString(frame_t);
						return ERROR_MSG;
					}
					ind = ind3D_to_ind(Xi, prm->Nslc);
					if(ind >= srf->frames[i].regRays.size()){
						error_handl_string = "ind = " + toString(ind) + "; rays.size = " + toString(srf->frames[i].regRays.size()) + "; Nslc = " + toStr(prm->Nslc);
						return ERROR_MSG;
					}
					srf->frames[i].regRays[ind].add(this);
					break;
				}
				frame_t += prm->dt;
				++i;
			}
		}

		this->t += dlt_t;
		if(prm->use_det){
			if(trng->detector){ // found triangle is a part of a detector
				trng->detector->regRays.push_back(RegisteredRay(this));
			} // TODO frames - building registered curve on the fly
			if(trng->is_absorber){
				rays->quit_ray();
				return 0;
			}
		}

		rx = this->r + dlt_r; // collision point found
		n = trng->n;
	}

	// phi - angle of incidence
	double cos_phi = fabs(dot(n,this->v)/(length(n)*length(this->v)));

	// ------------------------------ phys angles & amplitude ---------------------------
	double sin_phi = cos_sin(cos_phi);
	double p = sin_phi / this->c;
	//double p2 = p*p;
	double sin_p, sin_s;

	// snell's law
	// TODO sin > 1
	Material *mat_from = trng->mat[ bool2int(trng->sg(trng->r[0] - this->r, 0) < 0) ];
	sin_p = p * mat_from->Cp;
	sin_s = p * mat_from->Cs;

	// --------------------------------- post-geom - create new rays--------------------------------------
	double3 nr = normalize(cross(n,this->v)); // reflection surface
	double3 v_new, r_new;
	double A_new, t_new = this->t + prm->eps;

	double sh_abs = 0, sv_abs = 1;
	if(this->type == SRayType){
		sh_abs = dot(this->polar, nr);
		sv_abs = cos_sin(sh_abs);
	}

	// ---- build new P-ray -------
	A_new = this->A * sv_abs;
	if(sin_p <= 1){
		double cos_p = cos_sin(sin_p);
		v_new = newV(n, this->v, n*cos_p, cross(n, nr)*sin_p) * mat_from->Cp;
		if(is0(v_new)){
			error_handl_string = "v_new for P ray wasn't found for ray " + this->ToStr();
			return ERROR_MSG;
		}

		r_new = rx + v_new * prm->eps;

		rays->add_ray(new Ray(r_new, v_new, normalize(v_new), PRayType, A_new / sqrt(2), t_new));
	} else {
		rays->lostP += A_new * A_new;
	}
	// Ray(const double3 _r, const double3 _v, const int _type = BaseRayType, const double _A = 1, const double _t = 0):
	// newV(double3 n, double3 v, double3 vx, double3 vy)

	// ---- build new S-rays -------
	A_new = this->A;
	if(sin_s <= 1){
		double cos_s = cos_sin(sin_s);
		v_new = newV(n, this->v, n*cos_s, cross(n, nr)*sin_s) * mat_from->Cs;
		if(is0(v_new)){
			error_handl_string = "v_new for S ray wasn't found for ray " + this->ToStr();
			return ERROR_MSG;
		}

		if(sin_p > 1) // r_new wasn't assigned before
			r_new = rx + v_new * prm->eps;

		rays->add_ray(this->type == PRayType ?
					  // this->type == PRayType
					  new Ray(r_new, v_new, normalize(cross(v_new, nr)), SRayType, A_new, t_new) :
					  // this->type == SRayType
					  new Ray(r_new, v_new, normalize((this->polar + nr * sh_abs*(sqrt(2) - 1)) / sqrt(2)), SRayType,
							  A_new*sqrt(sh_abs*sh_abs + sv_abs*sv_abs/2), t_new));
		              // it can be done my creating 2 rays with || and _|_ polarisations, but it's equivalent to a single ray with a rotated polarisation);
	} else {
		rays->lostS += A_new * A_new;
	}

	rays->quit_ray();
	return 0;
}

pair<double, int> Ray::find_collision(Surface* srf, Params* prm)
{
	int i;
	double t, t_min = prm->Tmax * 2;
	Triangle *trngl;
	int i_coll = -1; // indicator of no found collision

	for(i = 0; i < srf->Npol; ++i){ // find collision point

		trngl = &(srf->polygons[i]); // so we don't have to call [i] every time. Also it's shorter
		t = dot(trngl->r[0] - this->r, trngl->n) / dot(this->v, trngl->n); // find time of collision

		if((SYS_EPS < t) && (t < t_min)){
		// if the possible collision can happen (t > 0) and if it's better than the one we already have (t < t_min)
			if(trngl->isInside(this->r + this->v*t, prm->eps)){ // if it's really the point, then save it
				t_min = t;
				i_coll = i;
			}
		}
	}

	return make_pair(t_min, i_coll);
}

string Ray::ToStr(const string spr)
{
	return spr + (this->type == PRayType ? "P-type" : "S-type") +
			"\nc = " + toString(this->c) + "; A = " + toString(this->A) + "; t = " + toString(this->t) +
			"\nr = " + toStr(this->r) +
			"\nv = " + toStr(this->v) + spr;
}

void Ray::print(ostream &output, string spr)
{
	output << this->ToStr(spr);
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ RayFront ----------------------------------------------
// --------------------------------------------------------------------------------------------

void RaysFront::add_ray(Ray *new_ray)
{
	this->last_ray->next = new_ray;
	this->last_ray = this->last_ray->next;
	this->inc_rays();
}

void RaysFront::shift_current_ray(void)
{
	Ray *old_ray = this->current_ray;
	this->current_ray = this->current_ray->next;
	delete old_ray;
}

void RaysFront::quit_ray(void)
{
	--this->n_alive_rays;
	this->shift_current_ray();
}

void RaysFront::inc_rays(void)
{
	++this->n_alive_rays;
	++this->n_total_rays;
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Surface -----------------------------------------------
// --------------------------------------------------------------------------------------------

Surface::Surface(int _n_pol, int _n_mat, int _n_det)
{
	this->resize_clear(_n_pol, _n_mat, _n_det);
}

void Surface::resize_clear(int _n_pol, int _n_mat, int _n_det)
{
	if(_n_pol != this->Npol){
		delete[] this->polygons;
		this->polygons = (_n_pol == 0 ? nullptr : (new Triangle[_n_pol]));
		this->Npol = _n_pol;
	}
	if(_n_mat != this->Nmat){
		delete[] this->materials;
		this->materials = (_n_mat == 0 ? nullptr : (new Material[_n_mat]));
		this->Nmat = _n_mat;
	}
	if(_n_det != this->Ndet){
		delete[] this->detectors;
		this->detectors = (_n_det == 0 ? nullptr : (new Detector[_n_det]));
		this->Ndet = _n_det;
	}
}

void Surface::print(ostream &output, string spr1, string spr2)
{
	int i;
	output << spr1 << "---------- Materials ----------\n"
		   << "Nmat = " << this->Nmat << "\n";
	for(i = 0; i < this->Nmat; ++i){
		this->materials[i].print(output, toString(i+1) + ") ");
	}
	output << "---------- Polygons ----------\n"
		   << "Npol = " << this->Npol << "\n";
	for(i = 0; i < this->Npol; ++i){
		this->polygons[i].print(output, "\n------------\n" + toString(i+1) + "\n");
	}
	output << spr2;
}

int Surface::load_from_file(string surf_filename, string mat_filename, Params *prm)
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
		input.close();
		error_handl_string = "Nmat = " + toString(_n_mat) + "; must be at least 2\n";
		return LESS_2_MATERIALS;
	}
	this->resize_clear(this->Npol, _n_mat, this->Ndet);
	for(i = 0; i < this->Nmat; ++i){
		input >> this->materials[i].Cp >> this->materials[i].Cs;
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

	// --------------------------------------- read vertices -----------------------------------------------
	double3 *vertex;
	int Nvertex;

	input >> Nvertex;
	vertex = new double3[Nvertex];
	bool find_Xbounds = areEq(prm->Xmax, prm->Xmin);
	if(find_Xbounds){
		prm->Xmax = make_double3(-9999999999999.0, -9999999999999.0, -9999999999999.0);
		prm->Xmin = make_double3(9999999999999.0, 9999999999999.0, 9999999999999.0);
	}

	for(i = 0; i < Nvertex; ++i){
		input >> vertex[i].x >> vertex[i].y >> vertex[i].z;

		if(find_Xbounds){ // find global borders
			if(vertex[i].x > prm->Xmax.x) prm->Xmax.x = vertex[i].x;
			if(vertex[i].x < prm->Xmin.x) prm->Xmin.x = vertex[i].x;
			if(vertex[i].y > prm->Xmax.y) prm->Xmax.y = vertex[i].y;
			if(vertex[i].y < prm->Xmin.y) prm->Xmin.y = vertex[i].y;
			if(vertex[i].z > prm->Xmax.z) prm->Xmax.z = vertex[i].z;
			if(vertex[i].z < prm->Xmin.z) prm->Xmin.z = vertex[i].z;
		}
	}

	// now we know real Xmax & Xmin & dX
	// so we can allocate memory for frames (and for the grid if necessary)
	if(prm->dX.x == 0) prm->dX.x = (prm->Xmax.x - prm->Xmin.x)*1.01; // bigger than max delta, so int(delta/dX) == 0 ,
	if(prm->dX.y == 0) prm->dX.y = (prm->Xmax.y - prm->Xmin.y)*1.01; // so Nslc = 1 in the end
	if(prm->dX.z == 0) prm->dX.z = (prm->Xmax.z - prm->Xmin.z)*1.01;
	prm->Nslc = get_Nslc(prm->Xmax - prm->Xmin, prm->dX) + 1;

	if(prm->draw_mov){
		this->frames = new Frame[prm->Nfrm + 1];
		if(prm->prnt_mode == FRAMES_DATA_MODE){
			int Nnodes = prm->Nslc.x * prm->Nslc.y * prm->Nslc.z;
			int ind;
			int3 Xi;

			for(i = 0; i < prm->Nfrm; ++i){
				this->frames[i].regRays.resize(Nnodes);
				/*
				 * now (27.07.2018) sizeof(Ray)==112, so here we need a lot of RAM
				 * In fact we don't need fields "c","t","next" for registered rays
				 * so we can create a special class for frame-registered rays which would be ~90 bytes.
				 * It's not a big difference, so I didn't bother so far, but it can be done at any moment.
				 */
				for(Xi.z = 0; Xi.z < prm->Nslc.z; ++Xi.z) for(Xi.y = 0; Xi.y < prm->Nslc.y; ++Xi.y) for(Xi.x = 0; Xi.x < prm->Nslc.x; ++Xi.x){
					ind = ind3D_to_ind(Xi, prm->Nslc);
					this->frames[i].regRays[ind].t = i * prm->dt;
					this->frames[i].regRays[ind].r = (Xi + 0.5) * prm->dX + prm->Xmin;
				}
			}
		}
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
			if(i_vertex[i2] >= Nvertex){
				error_handl_string = "index of the " + toString(i2) + "th vertex of the " + toString(i) + "th tetrahedron is " + toString(i_vertex[i2]) + "; Nvertex = " + toString(Nvertex) + "\n";
				delete[] vertex;
				read_polygons.clear();
				input.close();
				return ERROR_MSG;
			 }
			// TODO check for identical vertices in a single tetrahedron
		}
		input >> i_mat; // read material of tetrahedron_i
		// materils are indexed from 0, but 0th material is for background, so we don't do --i_mat
		if(i_mat >= this->Nmat){
			error_handl_string = "index of material of the " + toString(i) + "th tetrahedron is " + toString(i_mat) + "; Nmat = " + toString(this->Nmat) + "\n";
			delete[] vertex;
			read_polygons.clear();
			input.close();
			return ERROR_MSG;
		}

		for(i2 = 0; i2 < 4; ++i2){                        // for each tetrahedron sub-surface (tet_surf_ind[i2][:]) check overlaps
			for(i3 = 0; i3 < read_polygons.size(); ++i3){ // with all existing ones
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
				return ERROR_MSG;
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
			if(!read_polygons[i].mat[i2]){
				read_polygons[i].mat[i2] = this->materials; // &(material[0])
			}
		}
		if(read_polygons[i].mat[0]->isEq(read_polygons[i].mat[1])){
			// read_polygons[i] is a formal border with no real meaning
			// actually it can't happen here if everything before worked fine.
			read_polygons.erase(read_polygons.begin()+i);

			error_handl_string = toString(i) + "-th polygon has no material set\n";
			delete[] vertex;
			read_polygons.clear();
			input.close();
			return ERROR_MSG;
		}
	}

	this->resize_clear(read_polygons.size(), this->Nmat, this->Ndet);
	for(i = 0; i < this->Npol; ++i){ // copy read data
		this->polygons[i] = read_polygons[i];
	}
	read_polygons.clear();
	delete[] vertex;

	// --------------------------------------- read detectors -----------------------------------------------
	if(prm->use_det){
		int Ntrg, i_det, _n_det;
		Triangle *trng;

		input >> Ntrg >> _n_det;
		this->resize_clear(this->Npol, this->Nmat, _n_det);
		for(i = 0; i < Ntrg; ++i){
			for(i2 = 0; i2 < 3; ++i2){
				input >> i_vertex[i2];
				--i_vertex[i2];
				if(i_vertex[i2] >= Nvertex){
					error_handl_string = "index of the " + toString(i2) + "-th vertex of the " + toString(i) + "-th tetrahedron is " + toString(i_vertex[i2]) + "; Nvertex = " + toString(Nvertex) + "\n";
					input.close();
					return ERROR_MSG;
				 }
			}
			input >> i_det;
			--i_det;
			// --i is necessary because all arrays are indexed from 0 but all abjects in the file are indexed from 1
			if(i_det >= this->Ndet){
				error_handl_string = "index of the detector of the " + toString(i) + "-th polygon is " + toString(i_det) + "; Ndet = " + toString(this->Ndet) + "\n";
				input.close();
				return ERROR_MSG;
			}

			trng = this->findPolygon(i_vertex);
			if(trng){
				if(i_det >= 0){
					trng->detector = &(this->detectors[i_det]); // add i2-th polygon to the i_det-th detector
				}
				trng->is_absorber = 1; // just an absorbing polygon, not a detector
			} else {
				error_handl_string = toString(i) + "-th polygon is missing in the final version of the surface\n";
				input.close();
				return ERROR_MSG;
			}
		}
	}

	input.close();

	return 0;
}

Triangle* Surface::findPolygon(int* i_vertex) // find triangle by its vertices
{
	int i, i2, i3;
	int n_overlap_vrtx;

	for(i = 0; i < this->Npol; ++i){     // for each polygon
		n_overlap_vrtx = 0;
		for(i2 = 0; i2 < 3; ++i2){          // i_vertex[i2]
			for(i3 = 0; i3 < 3; ++i3){      // this->polygons[i].i_vrtx[i3]
				if(i_vertex[i2] == this->polygons[i].i_vrtx[i3]){
					++n_overlap_vrtx; // count overlapping vertices
					//cout << i2 << " " << i << " " << i3 << "\n";
				}
			}
		}
		if(n_overlap_vrtx >= 3){  // i - found polygon index
			if(n_overlap_vrtx > 3){
				error_handl_string = "Polygon (" + toString(i_vertex[0]) + ";" + toString(i_vertex[0]) + ";" + toString(i_vertex[0]) + ") has too many overlaps with " + toString(i) + "th polygon";
				CHECK(ERROR_MSG);
				exit(1);
			}
			break;
		}
	}

	// i == this->Npol means no existing polygon matches all 3 vertices
	return (i == this->Npol ? nullptr :  &(this->polygons[i]));
}

int Surface::saveDetectorInfo(string filename, Params *prm)
{
	ofstream output;

	output.open(filename);
	if(!output){
		output.close();
		error_handl_string = "can not open output file " + filename + " for writing\n";
		return CANT_OPEN_FILE_FOR_WRITING;
	}

	int i;
	// prm->Nfrm = (int)(prm->Tmax*(1 + prm->eps) / prm->dt)+1;
	output << this->Ndet << " " << prm->Nfrm << "\n";
	for(i = 0; i < this->Ndet; ++i){
		//output << i+1 << "\n";
		this->detectors[i].saveInfo(output, prm);
	}

	output.close();
	return 0;
}

int Surface::saveMovie(Params *prm)
{
	string path = "./" + prm->model_name + "/frames/";
    int i;
	ofstream output;
	string filename;

	time_t real_start_t;

	double tot_rays = 0, rays_printed = 0;
	for(i = 0; i < prm->Nfrm; ++i){
		tot_rays += this->frames[i].regRays.size();
	}

	time(&real_start_t);
	for(i = 0; i < prm->Nfrm; ++i){
		filename = path + toString(i) + ".frm";
		output.open(filename);
		if(!output){
			output.close();
			error_handl_string = "can not open output file \n|" + filename + "|\nfor writing\n";
			return CANT_OPEN_FILE_FOR_WRITING;
		}

		this->frames[i].saveToFile(output, prm);
		rays_printed += this->frames[i].regRays.size();

		output.close();

		time_progress(real_start_t, time(0), rays_printed / tot_rays, "saving movie");
	}

	return 0;
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Frame -------------------------------------------------
// --------------------------------------------------------------------------------------------

int Frame::saveToFile(ostream &output, Params *prm)
{
	int i;

	switch(prm->prnt_mode){ // TODO PRNT_MODE_ID
	case RAW_DATA_MODE:
	case FRAMES_DATA_MODE:
		int real_N = 0;
		for(i = 0; i < this->regRays.size(); ++i){
			if(this->regRays[i].A > 0)
				++real_N;
		}

		Ray *curr_ray;
		output << real_N << "\n";
		for(i = 0; i < this->regRays.size(); ++i){
			curr_ray = &(this->regRays[i]);
			if(curr_ray->A > 0){
				output //<< i << " "
				   	   << (curr_ray->type - BaseRayType) << " "
				   	   << curr_ray->t << " "
				   	   << curr_ray->A << " "
				   	   << curr_ray->c << " "
				   	   << toStr(curr_ray->r, "", " ", "") << " "
				   	   << toStr(curr_ray->v, "", " ", "") << "\n";
			}
		}
		break;
	}

	return 0;
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Detector ----------------------------------------------
// --------------------------------------------------------------------------------------------

void Detector::saveInfo(ostream &output, Params *prm)
{
	int i;
	switch(prm->prnt_mode){
		case RAW_DATA_MODE:
			// save raw data
			for(i = 0; i < this->regRays.size(); ++i){
				output << i << " " << this->regRays[i].t << " " << this->regRays[i].A << " " << this->regRays[i].c << " " << this->regRays[i].type-RayTimeIND << "\n";
			}
			break;
		case FRAMES_DATA_MODE:
			// save registered curve
			double *regArr = new double[prm->Nfrm];

#ifdef _OPENMP
#pragma omp parallel for
#endif
			for(i = 0; i < prm->Nfrm; ++i){
				regArr[i] = this->regValue(prm->dt*i, prm);
			}

			for(i = 0; i < prm->Nfrm; ++i){
				output << prm->dt*i << " " << regArr[i] << "\n";
			}

			delete[] regArr;
	}
}

double Detector::regValue(double t0, Params *prm)
{
	int i;
	double res = 0;

	for(i = 0; i < this->regRays.size(); ++i){
		// res += this->peakFnc((t0 - this->regRays[i].t) / prm->tau, prm->B, this->regRays[i].A);
		res += this->peakFnc((t0 - this->regRays[i].t) / prm->tau, prm->B);
	}

	return res / prm->Nrays;
}

double Detector::peakFnc(double x, double B, double A)
// y = A*cos(2pi*f*t)*cos^4(t/tau) =
// = A*cos(B*x)*cos^4(x), x = (t-t0)/tau, B = f*tau
{
	x = std::abs(x);
	if(x > pi_d2){
		return 0;
	} else {
		double _b = float_part(B*x);
		if((0.25 <= _b) && (_b <= 0.75)) // cos(B*x) < 0
			return 0;
		_b = cos(x);
		_b = _b * _b; // cos^2
		_b = A * cos(pi_m2 * B * x) * _b * _b;
		return  _b;
	}
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Params ------------------------------------------------
// --------------------------------------------------------------------------------------------

int Params::load_from_file(string filename)
{
	ifstream input(filename);
	if(!input){
		input.close();
		error_handl_string = "file " + filename + " is missing\n";
		return ERROR_MSG;
	}

	string buf_s;
	std::getline(input, buf_s); // read comment line
	input >> this->Nrays >> this->Tmax >> this->dt >> this->Amin >> this->f >> this->tau >> this->eps
	      >> this->use_det >> this->draw_mov >> this->prnt_mode
	      >> this->Xmin.x >> this->Xmin.y >> this->Xmin.z
	      >> this->Xmax.x >> this->Xmax.y >> this->Xmax.z
	      >> this->dX.x >> this->dX.y >> this->dX.z;
	this->prnt_mode += PRINT_MODE_ID;

	input.close();

	if(this->eps == 0) this->eps = SYS_EPS;
	if(this->dt > 0){
		this->Nfrm = (int)(this->Tmax*(1 + this->eps) / this->dt)+1;
	} else {
		this->Nfrm = 0;
		if((this->prnt_mode == FRAMES_DATA_MODE) && this->use_det){
			error_handl_string = "average-sum result for detectors is requested, but dt <= 0\n";
			return ERROR_MSG;
		}
		if(draw_mov){
			error_handl_string = "draw movie requested, but dt <= 0\n";
			return ERROR_MSG;
		}
	}

	this->B = this->f * this->tau;
	this->f *= TIME_UNIT;   //  Gz * 10^6
	this->tau /= TIME_UNIT; // sec * 10^-6

	return 0;
}

void Params::print_full(ostream &output, string spr)
{
	for(int i = 0; i < this->paramFHead.size(); ++i){
		output << this->paramFHead[i] << " | ";
	}
	output << spr << this->Nrays << spr
		   << this->Tmax << spr
		   << this->Amin << spr
		   << this->eps << spr
		   << "model_name: |" << this->model_name << "|" << spr;
		   //<< "alive_rays: " << this->n_alive_rays << "; total_rays: " << this->n_total_rays << spr;
}

void Params::print(ostream &output)
{
	int spForVal = 15;

	for(int i = 0; i < this->paramFHead.size(); ++i){
		output << setw(spForVal) << this->paramFHead[i];
	}

	output << "\n"
		   << setw(spForVal) << this->Nrays << setw(spForVal) << this->Tmax << setw(spForVal) << this->dt
		   << setw(spForVal) << this->Amin << setw(spForVal) << this->f/TIME_UNIT << setw(spForVal) << this->tau*TIME_UNIT
		   << setw(spForVal) << this->eps << setw(spForVal) << this->use_det << setw(spForVal) << this->draw_mov
		   << setw(spForVal) << (this->prnt_mode - PRINT_MODE_ID)
		   << setw(spForVal) << this->Xmin.x << setw(spForVal) << this->Xmin.y << setw(spForVal) << this->Xmin.z
		   << setw(spForVal) << this->Xmax.x << setw(spForVal) << this->Xmax.y << setw(spForVal) << this->Xmax.z
		   << setw(spForVal) << this->dX.x << setw(spForVal) << this->dX.y << setw(spForVal) << this->dX.z;
}

int Params::save_to_file(string filename)
{
	ofstream output(filename);
	if(!output){
		error_handl_string = "can not open file                                    \n|" +
				              filename +
				              "|                                             \nfor writing                                           \n";
		return CANT_OPEN_FILE_FOR_WRITING;
	}

	this->print(output);

	output.close();
	return 0;
}

// --------------------------------------------------------------------------------------------
// ------------------------------------ Global Fncs -------------------------------------------
// --------------------------------------------------------------------------------------------

int compute(Surface *srf, Params *prm, RaysFront *rays)
{
	time_t real_start_t;
	int i;

	double3 *rays_v = new double3[prm->Nrays];
	// generate rays
	for(i = 0; i < prm->Nrays; ++i){
		rays_v[i] = vecByAngles(myRnd(pi/4, 3*pi/4), 0) * srf->materials[1].Cp; // v || oY
		//rays_v[i] = vecByAngles(myRnd(pi/4, 3*pi/4), myRnd(-pi/5, pi/5)) * srf->materials[1].Cp; // v || oY
	}
	/*
	for(i = 0; i < prm->Nrays/4; ++i){
		// rays_v[4*i] = rndVec(srf->materials[1].Cp);
		rays_v[4*i]   = vecByAngles(myRnd(-pi/4-0.1, pi/4+0.1), myRnd(-pi/4-0.1, pi/4+0.1)) * srf->materials[1].Cp;
		rays_v[4*i+1] = vecByAngles(myRnd(pi/4-0.1, pi*3/4+0.1), myRnd(-pi/4-0.1, pi/4+0.1)) * srf->materials[1].Cp;
		rays_v[4*i+2] = vecByAngles(myRnd(pi*3/4-0.1, pi*5/4+0.1), myRnd(-pi/4-0.1, pi/4+0.1)) * srf->materials[1].Cp;
		rays_v[4*i+3] = vecByAngles(myRnd(pi*5/4-0.1, pi*7/4+0.1), myRnd(-pi/4-0.1, pi/4+0.1)) * srf->materials[1].Cp;
		// vecByAngles(myRnd(-pi, pi), myRnd(-pi_d2, pi_d2))*V
	}
	*/

	time(&real_start_t);
	for(i = 0; i < prm->Nrays; ++i){
		rays->current_ray = new Ray(make_double3(0,0,0), rays_v[i], normalize(rays_v[i]), PRayType, 1, 0);
		// Ray(const double3 _r, const double3 _v, const double3 _polar, const int _type = PRayType, const double _A = 1, const double _t = 0, const Ray *_next = nullptr);
		rays->inc_rays();
		rays->last_ray = rays->current_ray;
		rays->totalE0 += rays->current_ray->A * rays->current_ray->A;

		do{
			if(rays->current_ray->move(srf, prm, rays)){
				error_handl_string += ("\n" + toString(i) + "-th ray failed\n");
				delete[] rays_v;
				return ERROR_MSG;
			}
		}while(rays->current_ray); // while there are rays to compute

		time_progress(real_start_t, time(0), (i+1) / (float(prm->Nrays)), "computing");
	}

	cout << "lostP/totalE = " << rays->lostP / rays->totalE0 << "                                      \n"
		 << "lostS/totalE = " << rays->lostS / rays->totalE0 << "                                      \n";

	delete[] rays_v;
	return 0;
}
