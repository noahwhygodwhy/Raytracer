#include "shape.cl"
#include "sharedStructs.cl"

HitResult rayHitTriangle(const Triangle tri, const Ray ray, __global const Vertex* vertices, int shapeIdx) {

	
    HitResult hit;
    hit.hit = false;


	//TODO: should really be done every frame, not every frame*ray

	// constexpr float epsilon = epsilon<float>();
	float3 a = vertices[tri.vertA].position.xyz;
	float3 b = vertices[tri.vertB].position.xyz;
	float3 c = vertices[tri.vertC].position.xyz;

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

	float3 badNormal = ((vertices[tri.vertA].normal * bary.x) + (vertices[tri.vertB].normal * bary.y) + (vertices[tri.vertC].normal * bary.z)).xyz;


	hit.normal = badNormal;// transformNormal(badNormal, model);



	hit.uv = ((vertices[tri.vertA].uv * bary.x) + (vertices[tri.vertB].uv * bary.y) + (vertices[tri.vertC].uv * bary.z)).xy;
	
	hit.matIdx = tri.shape.matIdx;


	//barycentric.z = dot(tri.E2, qvec);

	return hit;
}