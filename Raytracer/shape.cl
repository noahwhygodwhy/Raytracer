#ifndef SHAPE_H
#define SHAPE_H

#include "sharedStructs.cl"
#include "helpers.cl"



enum Side {
	RIGHT,
	LEFT,
	MIDDLE
};



rayAABBResult rayAABB(AABB this, const Ray ray) { 

    rayAABBResult result;
	float3 invD = 1.0f / ray.direction.xyz;
	float3 t0 = (this.min.xyz - ray.origin.xyz) * invD;
	float3 t1 = (this.max.xyz - ray.origin.xyz) * invD;

	float3 tSmaller = min(t0, t1);
	float3 tBigger = max(t0, t1);

	float tMin = -INFINITY;
	float tMax = INFINITY;
	tMin = max(tMin, max(tSmaller.x, max(tSmaller.y, tSmaller.z)));
	tMax = min(tMax, min(tBigger.x, min(tBigger.y, tBigger.z)));

	result.hit = tMin < tMax && tMax >= 0.0;
	result.enter = ray.origin + (ray.direction * tMin);
	result.exit = ray.origin + (ray.direction * tMax);

	return result;
}

#endif