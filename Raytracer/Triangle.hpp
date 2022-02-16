#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Ray.hpp"
#include "Shape.hpp"
#include "Mesh.hpp"

//extern uint32_t frameX;
//extern uint32_t frameY;

extern bool prd;



float antiNDC(float i) {
	return (i + 1.0f) / 2.0f;
}



float triArea(vec3 a, vec3 b, vec3 c) {
	vec3 d = b - a;
	vec3 e = c - a;
	vec3 f = glm::cross(d, e);
	return length(f) / 2.0f;
}

class Triangle : public Shape {
	Vertex vertices[3];

	//uvec2 zLimits;
	bool onscreen = true;
	//dvec3 E1;
	//dvec3 E2; 
	//float planeD;

public:
	//dvec3 maxBounding;
	//dvec3 minBounding;
	//dvec3 middleBounding;
	//dvec3 normal;
	Mesh* mesh;
	//Triangle() {}

	Triangle(Vertex a, Vertex b, Vertex c, Mesh* m)
		: Shape(AABB(glm::max(a.position, glm::max(b.position, c.position)), 
			glm::min(a.position, glm::min(b.position, c.position))), materials.at("Copper"))
	{
		this->mesh = m;
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;

		//this->middleBounding = (this->maxBounding + this->minBounding) / 2.0;

		//E1 = b.position - a.position;
		//E2 = c.position - a.position;

		/*this->normal.x = (E1.y * E2.z) - (E1.z * E2.y);
		this->normal.y = (E1.z * E2.x) - (E1.x * E2.z);
		this->normal.z = (E1.x * E2.y) - (E1.y * E2.x);
		this->normal = glm::normalize(this->normal);*/


		//this->planeD = glm::dot(this->normal, this->vertices[0].position);
	}

	vec3 getP(const mat3x3& transformedVerts, vec3 bary) const {
		return mat3x3(vertices[0].position, vertices[1].position, vertices[2].position) * bary;
	}
	Vertex operator[](uint32_t index) const {
		return vertices[index];
	}
	bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const;
	dvec3 getColor(const HitResult& hit)const;

private:

};


bool Triangle::rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime) const {


	//TODO: should really be done every frame, not every frame*ray
	constexpr double epsilon = glm::epsilon<double>();
	dvec3 a = transformPos(this->vertices[0].position, dmat4(1.0), view);
	dvec3 b = transformPos(this->vertices[1].position, dmat4(1.0), view);
	dvec3 c = transformPos(this->vertices[2].position, dmat4(1.0), view);

	dvec3 E1 = b - a;
	dvec3 E2 = c - a;

	dvec3 transformedNormal;

	transformedNormal = transformNormal(glm::cross(E1, E2), dmat4(1));

	dvec3 pvec = glm::cross(ray.direction, E2);
	double determinant = glm::dot(E1, pvec);

	if (determinant < epsilon) {//TODO: i fucked something up
		return false;
	}

	double inverseDeterminant = 1.0 / determinant;

	dvec3 tvec = ray.origin - a;

	double u = glm::dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0 || u > 1.0) {
		return false;
	}

	dvec3 qvec = glm::cross(tvec, E1);
	double v = glm::dot(ray.direction, qvec) * inverseDeterminant;

	if (v < 0.0 || u + v > 1.0) {
		return false;
	}


	hit.bary = dvec3(1.0f - (u + v), u, v);
	hit.position = mat3x3(a, b, c) * hit.bary;

	if (prd)printf("hit at %s\n",glm::to_string(hit.position).c_str());

	hit.depth = glm::length(hit.position - ray.origin);
	hit.normal = glm::normalize(transformedNormal);
	hit.uv = mat3x2(this->vertices[0].texCoords, this->vertices[1].texCoords, this->vertices[2].texCoords)* hit.bary;
	hit.shape = (Triangle*)this;



	//barycentric.z = glm::dot(tri.E2, qvec);

	return true;
}
//TODO: make this use the material
dvec3 Triangle::getColor(const HitResult& hit)const {
	dvec2 squares(20.0, 40.0);
	//double squares = 10.0;

	ivec2 flatUV = glm::floor(hit.uv * squares);
	//int32_t flatU = int32_t(floor(hit.uv.x * squares));
	//int32_t flatV = int32_t(floor(hit.uv.y * squares));
	//printf("hit.uv: %s\n", glm::to_string(hit.uv).c_str());
	//printf("flatU: %u, flatV: %u\n", flatU, flatV);
	//printf("%s\n", (flatU + flatV) % 2 == 0 ? "true" : "false");

	//return dvec3(hit.uv.x, hit.uv.y, 0.0);

	if ((flatUV.x + flatUV.y) % 2 == 0) {
		//printf("color 1\n");
		return dvec3(1.0, 0.0, 0.0);

	}
	else {
		//printf("color 2======================\n");
		return dvec3(1.0, 1.0, 0.0);
	}
}

float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}
//sign() and pointInTriangle() based off of the code in this stackoverflow post https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float sign(vec2 a, vec2 b, vec2 c) {
	return ((a.x - c.x) * (b.y - c.y)) - ((b.x - c.x) * (a.y - c.y));
}

/*enum TSide {
	TRIGHT,
	TLEFT,
	TMIDDLE
};
#define DIMS 3*/
//based on the fast ray-box intersection chapter by andrew woo from Graphics Gems
//calculates hits between a ray and a triangle's bounding box
/*bool rayHit(const Ray& ray, const vec3& minBounding, const vec3& maxBounding) {

	TSide sides[DIMS]; //quadrants
	double planes[DIMS]; //candiatePlane
	double maxT[DIMS];
	bool inside = true;

	for (uint32_t i = 0; i < DIMS; i++) {
		if (ray.origin[i] < minBounding[i]) {
			sides[i] = TLEFT;
			planes[i] = minBounding[i];
			inside = false;
		}
		else if (ray.origin[i] > maxBounding[i]) {
			sides[i] = TRIGHT;
			planes[i] = maxBounding[i];
			inside = false;
		}
		else {
			sides[i] = TMIDDLE;
		}
	}
	if (inside) {
		return true;
	}
	for (uint32_t i = 0; i < DIMS; i++)
	{
		if (sides[i] != TMIDDLE && ray.direction[i] != 0.) {
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
	vec3 coord;

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
}*/


//actually calculates if a ray hits the triangle, and returns the barycentric coordinate
/*bool rayHit(const Ray& ray, const Triangle& tri, dvec3& barycentric, dvec3& hitPoint) {

	dvec3 pvec = glm::cross(ray.direction, tri.E2);
	double determinant = glm::dot(tri.E1, pvec);


	if(determinant < glm::dot(tri.E1, pvec)) {
		return false;
	}

	double inverseDeterminant = 1.0f / determinant;

	dvec3 tvec = ray.origin - tri[0].position;

	double u = glm::dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0f || u > 1.0f) {
		return false;
	}

	dvec3 qvec = glm::cross(tvec, tri.E1);
	double v = glm::dot(ray.direction, qvec)* inverseDeterminant;

	if (v < 0.0f || u + v > 1.0f) {
		return false;
	}

	barycentric = vec3(1.0f - (u + v), u, v);
	hitPoint = tri.getP(barycentric);
	//barycentric.z = glm::dot(tri.E2, qvec);

	return true;
}*/


//actually calculates if a ray hits the triangle, and returns the barycentric coordinate
/*double kdRayHit(const Ray& ray, const Triangle& tri, vec3& barycentric) {

	dvec3 pvec = glm::cross(ray.direction, tri.E2);
	double determinant = glm::dot(tri.E1, pvec);

	if (determinant < glm::dot(tri.E1, pvec)) {
		return INFINITY;
	}

	double inverseDeterminant = 1.0f / determinant;

	dvec3 tvec = ray.origin - tri[0].position;

	double u = glm::dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0f || u > 1.0f) {
		return INFINITY;
	}

	dvec3 qvec = glm::cross(tvec, tri.E1);
	double v = glm::dot(ray.direction, qvec) * inverseDeterminant;

	if (v < 0.0f || u + v > 1.0f) {
		return INFINITY;
	}
	barycentric = vec3(1.0f - (u + v), u, v);

	dvec3 hitPos = tri.getP(barycentric);

	double distance = glm::length2(hitPos - ray.origin);

	//barycentric.z = glm::dot(tri.E2, qvec);

	return distance;
}*/

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