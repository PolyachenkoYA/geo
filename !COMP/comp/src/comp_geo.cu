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
		prm.res_filename = prm.model_name + ".res";
		global_logFname = prm.model_name + ".log";
		SAY_LOG("session " + sessionID + "\n");
		// -------------------------------------------
		// ------------ read parameters --------------
		// -------------------------------------------
		if(CHECK(prm.load_from_file(prm.prm_filename))) return 1;

		// -----------------------------------------------
		// ---------------- read surface -----------------
		// -----------------------------------------------
		Surface srf;
		if(CHECK(srf.load_from_file(prm.tet_filename, prm.material_filename, &prm))) return 1;

		// ----------------------------------------------------
		// --------------- conduct computation ----------------
		// ----------------------------------------------------
		RaysFront rays;
		if(CHECK(compute(&srf, &prm, &rays))) return 1;

		// --------------------------------------------
		// ------------ save results ------------------
		// --------------------------------------------
		string path = "./" + prm.model_name + "/frames/";
	    mkdir(prm.model_name.c_str(),S_IRWXU | S_IRWXG);
	    mkdir(path.c_str(),S_IRWXU | S_IRWXG);

		if(CHECK(prm.save_to_file("./" + prm.model_name + "/" + prm.prm_filename))) return 1;
		if(prm.use_det){
			if(CHECK(srf.saveDetectorInfo(prm.res_filename, &prm))) return 1;
		}
		if(prm.draw_mov){
			if(CHECK(srf.saveMovie(&prm))) return 1;
		}
	} else {
		cout << "Usage:\n./comp_geo      model_name\n";
	}

	return 0;
}

