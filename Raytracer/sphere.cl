#include "shape.cl"


#include "sharedStructs.cl"

HitResult rayHitSphere(const Sphere sphere, const Ray ray, uint shapeIdx){

    HitResult hit;
    hit.hit = false;

	float3 L = ray.origin.xyz - sphere.origin.xyz;
	float a = dot(ray.direction, ray.direction);
	float b = 2.0f * dot(ray.direction.xyz, L);
	float c = dot(L, L) - (sphere.radius*sphere.radius);
	float t0, t1;
	
    float2 temp = solveQuadratic(a, b, c);
    t0 = temp.x;
    t1 =  temp.y;

	if (isnan(t0) && isnan(t1)) {
		return hit;
	}

	if (t0 > t1) {
        float s = t0;
        t0 = t1;
        t1 = s;
	}
	if (t0 < 0.0f) {
		t0 = t1;
		if (t0 < 0.0f) {
		    return hit;
		}
	}
	
	hit.position = ray.origin + (ray.direction * t0);
	hit.normal = normalize(hit.position - sphere.origin);

	hit.uv = (float4)(atan2(hit.normal.x, hit.normal.z) / (2.0f * M_PI_F) + 0.5f, hit.normal.y * 0.5f + 0.5f, 0.0, 0.0);

	hit.depth = t0;
	hit.shapeIdx = shapeIdx;

	return hit;
}


//FLT_EPSILON     