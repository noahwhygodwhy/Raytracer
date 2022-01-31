#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Ray.hpp"


extern uint32_t frameX;
extern uint32_t frameY;



float antiNDC(float i) {
	return (i + 1.0f) / 2.0f;
}



float triArea(vec3 a, vec3 b, vec3 c) {
	vec3 d = b - a;
	vec3 e = c - a;
	vec3 f = glm::cross(d, e);
	return length(f) / 2.0f;
}

class Triangle {
	Vertex vertices[3];

	//uvec2 zLimits;
	bool onscreen = true;
	fvec3 U;
	fvec3 V; 
	//float planeD;

public:
	fvec3 maxBounding;
	fvec3 minBounding;
	fvec3 E1;
	fvec3 E2;
	vec3 normal;
	Mesh* mesh;
	Triangle() {}

	Triangle(Vertex a, Vertex b, Vertex c, Mesh* m) {
		this->mesh = m;
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;

		this->maxBounding = glm::max(a.position, glm::max(b.position, c.position));
		this->minBounding = glm::min(a.position, glm::min(b.position, c.position));

		U = b.position - a.position;
		V = c.position - a.position;

		this->normal.x = (U.y * V.z) - (U.z * V.y);
		this->normal.y = (U.z * V.x) - (U.x * V.z);
		this->normal.z = (U.x * V.y) - (U.y * V.x);
		this->normal = glm::normalize(this->normal);


		this->E1 = b.position - a.position;
		this->E2 = c.position - a.position;

		//this->planeD = glm::dot(this->normal, this->vertices[0].position);
	}
	Vertex operator[](uint32_t index) const {
		return vertices[index];
	}
private:
};
float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}
//sign() and pointInTriangle() based off of the code in this stackoverflow post https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float sign(vec2 a, vec2 b, vec2 c) {
	return ((a.x - c.x) * (b.y - c.y)) - ((b.x - c.x) * (a.y - c.y));
}

enum Side {
	RIGHT,
	LEFT,
	MIDDLE
};
#define DIMS 3
//based on the fast ray-box intersection chapter by andrew woo from Graphics Gems
bool rayHit(const Ray& ray, const vec3& minBounding, const vec3& maxBounding, vec3& coord) {

	Side sides[DIMS]; //quadrants
	float planes[DIMS]; //candiatePlane
	float maxT[DIMS];
	bool inside = true;

	for (uint32_t i = 0; i < DIMS; i++) {
		if (ray.origin[i] < minBounding[i]) {
			sides[i] = LEFT;
			planes[i] = minBounding[i];
			inside = false;
		}
		else if (ray.origin[i] > maxBounding[i]) {
			sides[i] = RIGHT;
			planes[i] = maxBounding[i];
			inside = false;
		}
		else {
			sides[i] = MIDDLE;
		}
	}
	if (inside) {
		return true;
	}
	for (uint32_t i = 0; i < DIMS; i++)
	{
		if (sides[i] != MIDDLE && ray.direction[i] != 0.) {
			maxT[i] = (planes[i] - ray.origin[i]) / ray.direction[i];
		}
		else {
			maxT[i] = -1.;
		}
	}

	uint32_t intersectionDim = 0;
	for (uint32_t i = 1; i < DIMS; i++) {
		if (maxT[intersectionDim] < maxT[i])
		{
			intersectionDim = i;
		}
	}

	if (maxT[intersectionDim] < 0.0f) { 
		return (false); 
	}

	for (uint32_t i = 0; i < DIMS; i++)
		if (intersectionDim != i) {
			coord[i] = ray.origin[i] + maxT[intersectionDim] * ray.direction[i];
			if (coord[i] < minBounding[i] || coord[i] > maxBounding[i])
				return false;
		}
		else {
			coord[i] = planes[i];
		}
	return true;
}


bool rayHit(const Ray& ray, const Triangle& tri, vec3& barycentric, vec3& hitPoint) {

	fvec3 pvec = glm::cross(ray.direction, tri.E2);
	float determinant = glm::dot(tri.E1, pvec);


	if(determinant < glm::dot(tri.E1, pvec)) {
		return false;
	}

	float inverseDeterminant = 1.0f / determinant;

	fvec3 tvec = ray.origin - tri[0].position;

	float u = glm::dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0f || u > 1.0f) {
		return false;
	}

	fvec3 qvec = glm::cross(tvec, tri.E1);
	float v = glm::dot(ray.direction, qvec)* inverseDeterminant;

	if (v < 0.0f || u + v > 1.0f) {
		return false;
	}

	barycentric = vec3(1.0f - (u + v), u, v);

	//barycentric.z = glm::dot(tri.E2, qvec);



	return true;
}

bool pointInTriangle(vec2 point, Triangle tri) {
	float signab, signbc, signca;
	bool hasNeg, hasPos;

	signab = sign(point, tri[0].position, tri[1].position);
	signbc = sign(point, tri[1].position, tri[2].position);
	signca = sign(point, tri[2].position, tri[0].position);

	hasNeg = (signab < 0) || (signbc < 0) || (signca < 0);
	hasPos = (signab > 0) || (signbc > 0) || (signca > 0);

	return !(hasNeg && hasPos);
}


#endif