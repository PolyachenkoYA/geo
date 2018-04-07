/*
 ============================================================================
 Name        : comp_geo.cu
 Author      : PolyachYA
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include "../../geo.cu"

int main(int argc, char **argv)
{

	if(argc >= 2){
		time_t now = time(0);
		string sessionID = string(ctime(&now));
		SAY_LOG("session " + sessionID + "\n");
		// ------------ read params --------------
		string params_filename = string(argv[1]);
		// prm is global. It's not good, but I don't see why.
		// So at least for now it's done this way in case that I don't know how to do better
		if(CHECK(prm.load_from_file(params_filename), params_filename)) return 1;
		if(argc > 2){
			prm.msh_filename = string(argv[2]);
		}

		// ------------ create test surface --------------
		Surface srf;
		if(CHECK(srf.load_from_file(prm.msh_filename), prm.msh_filename)) return 1;

		// ------------ create test rays ---------------
		//Wavefront waves;
		Ray ray(make_double3(0,0,0), normalize(make_double3(1,2,0)));
		//waves.add(ray);

		// ------------- conduct test computation -------
		ray.move(srf);
	} else {
		cout << "Usage:\n./comp_geo      params_file      [msh_file]\n";
	}

	return 0;
}

