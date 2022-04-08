#include "shape.cl"
#include "sharedStructs.cl"

HitResult rayHitTriangle(const UShape tri, const Ray ray, __global const Vertex* vertices, int shapeIdx) {

	
    HitResult hit;
    hit.hit = false;


	//TODO: should really be done every frame, not every frame*ray

	// constexpr float epsilon = epsilon<float>();
	Vertex vA = vertices[convert_int(tri.values.x)];
	Vertex vB = vertices[convert_int(tri.values.y)];
	Vertex vC = vertices[convert_int(tri.values.z)];


	float3 a = vA.position.xyz;
	float3 b = vB.position.xyz;
	float3 c = vC.position.xyz;

	float3 E1 = b - a;
	float3 E2 = c - a;

	float3 pvec = cross(ray.direction, E2);
	float determinant = dot(E1, pvec);

	if (determinant < 0.0f) {//TODO: i fucked something up
		return hit;
	}

	float inverseDeterminant = 1.0f / determinant;

	float3 tvec = ray.origin - a;

	float u = dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0f || u > 1.0f) {
		return hit;
	}

	float3 qvec = cross(tvec, E1);
	float v = dot(ray.direction, qvec) * inverseDeterminant;

	if (v < 0.0f || u + v > 1.0f) {
		return hit;
	}


	float3 bary = (float3)(1.0f - (u + v), u, v);

    hit.hit = true;
	hit.position = (a * bary.x) + (b * bary.y) + (c * bary.z);


	hit.depth = distance(hit.position, ray.origin);

	float3 badNormal = ((vA.normal * bary.x) + (vB.normal * bary.y) + (vC.normal * bary.z)).xyz;


	hit.normal = badNormal;// transformNormal(badNormal, model);



	hit.uv = ((vA.uv * bary.x) + (vB.uv * bary.y) + (vC.uv * bary.z)).xy;
	
	hit.matIdx = tri.matIdx;


	//barycentric.z = dot(tri.E2, qvec);

	return hit;
}