#include "spec_math.cuh"

int rightV(double3* v, double3 n, double3 v0)
{
	int i, res = -1;
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

	return res;
}

double3 newV(double3 n, double3 v_old, double3 vx, double3 vy)
{
	double3 v_news[4];
	// we don't know the direction of n and nr, so we have to check all 4 possible variants
	// to choose the one that really is the physicaly reflected ray
	v_news[0] = vx + vy;
	v_news[1] = vx - vy;
	v_news[2] = -v_news[1];
	v_news[3] = -v_news[0];

	int i = rightV(v_news, n, v_old);
	return i == -1 ? make_double3(0,0,0) : v_news[i];
}
