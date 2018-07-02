/*
 ============================================================================
 Name        : comp_geo.cu
 Author      : PolyachYA
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include "../../geo.cuh"

int main(int argc, char **argv)
{

	if(argc >= 2){
		// --------------------------------------------
		// ------------ handle input ------------------
		// --------------------------------------------
		time_t now = time(0);
		string sessionID = string(ctime(&now));
		Params prm;
		prm.model_name = string(argv[1]);
		prm.prm_filename = prm.model_name + ".prm";
		prm.tet_filename = prm.model_name + ".tet";
		prm.material_filename = prm.model_name + ".mat";
		global_logFname = prm.model_name + ".log";
		SAY_LOG("session " + sessionID + "\n");
		// -------------------------------------------
		// ------------ read parameters --------------
		// -------------------------------------------
		if(CHECK(prm.load_from_file(prm.prm_filename), prm.prm_filename)) return 1;

		// -----------------------------------------------
		// ------------ create test surface --------------
		// -----------------------------------------------
		Surface srf;
		if(CHECK(srf.load_from_file(prm.tet_filename, prm.material_filename), prm.tet_filename + " or " + prm.material_filename)) return 1;

		// ---------------------------------------------
		// ------------ create test rays ---------------
		// ---------------------------------------------
		Ray ray(make_double3(0,0,0), normalize(make_double3(1,0,2)) * (srf.materials[1].Cp), PRayType);
		// Ray(const double3 _r, const double3 _v, const int _type = BaseRayType, const double _A = 1, const double _t = 0):

		// ----------------------------------------------
		// ------------- conduct test computation -------
		// ----------------------------------------------
		//prm.print(cout);
		//srf.print(cout);

		ray.move(&srf, &prm);
	} else {
		cout << "Usage:\n./comp_geo      model_name\n";
	}

	return 0;
}

