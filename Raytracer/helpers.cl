#ifndef HELPERS_H
#define HELPERS_H
#include "sharedStructs.cl"

float2 solveQuadratic(float a, float b, float c) {
	float discriminant = (b * b) - (4.0 * a * c);
	if (discriminant < 0.0) {
        return (float2)(NAN, NAN);
	}
	return (float2)( (-b + sqrt(discriminant)) / (2.0f * a), (-b - sqrt(discriminant)) / (2.0f * a));
}

float3 rotateVector(float3 v, float theta, float3 k) {
    return (v*cos(theta))+(cross(k, v)*sin(theta))+(k*dot(k, v)*(1.0f-cos(theta)));
}


#endif