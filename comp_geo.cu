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

int main(void)
{
	double h = 10;
	double big_dist = 2*MapSize;

	// ------------ create test surface --------------
	Surf surf;
	Triangle trng;
	double3 r1,r2,r3;
	// upper bound
	r1 = make_double3(-big_dist, h, big_dist);
	r2 = make_double3(-big_dist, h, -big_dist);
	r3 = make_double3(big_dist, h, big_dist);
	trng = Triangle(r1,r2,r3);
	surf.add(trng);
	r1 = make_double3(-big_dist, h, -big_dist);
	r2 = make_double3(big_dist, h, -big_dist);
	r3 = make_double3(big_dist, h, big_dist);
	trng = Triangle(r1,r2,r3);
	surf.add(trng);

	// lower bound
	r1 = make_double3(-big_dist, -h, big_dist);
	r2 = make_double3(-big_dist, -h, -big_dist);
	r3 = make_double3(big_dist, -h, big_dist);
	trng = Triangle(r1,r2,r3);
	surf.add(trng);
	r1 = make_double3(-big_dist, -h, -big_dist);
	r2 = make_double3(big_dist, -h, -big_dist);
	r3 = make_double3(big_dist, -h, big_dist);
	trng = Triangle(r1,r2,r3);
	surf.add(trng);

	// ------------ create test rays ---------------
	Wavefront waves;
	Ray ray(make_double3(0,0,0), normalize(make_double3(1,2,0)));
	waves.add(ray);

	// ------------- conduct test computation -------
	//move_ray(ray, surf);
	ray.move(surf);

	return 0;
}

