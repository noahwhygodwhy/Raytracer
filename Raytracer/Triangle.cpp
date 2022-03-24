#include "Triangle.hpp"


using namespace std;
using namespace glm;


bool Triangle::rayHit(const Ray& ray, HitResult& hit, double currentTime) const {

	

	if(prd)printf("ray checking triangle\n");


	//TODO: should really be done every frame, not every frame*ray

	constexpr double epsilon = glm::epsilon<double>();
	dvec3 a = this->vertices[0].position;
	dvec3 b = this->vertices[1].position;
	dvec3 c = this->vertices[2].position;


	if (prd)printf("a: %s\nb: %s\nc: %s\n", glm::to_string(a).c_str(), glm::to_string(b).c_str(), glm::to_string(c).c_str());


	if (prd)printf("a: %s\nb: %s\nc: %s\n", glm::to_string(this->vertices[0].normal).c_str(), glm::to_string(this->vertices[1].normal).c_str(), glm::to_string(this->vertices[2].normal).c_str());

	dvec3 E1 = b - a;
	dvec3 E2 = c - a;
	if (prd)printf("E1: %s\nE2: %s\n", glm::to_string(E1).c_str(), glm::to_string(E2).c_str());

	if (prd)printf("ray.direction: %s\n", glm::to_string(ray.direction).c_str());

	dvec3 pvec = glm::cross(ray.direction, E2);
	if (prd)printf("pvec: %s\n", glm::to_string(pvec).c_str());
	double determinant = glm::dot(E1, pvec);

	if (determinant < 0.0) {//TODO: i fucked something up
		if (prd) printf("triangle no hit 1\n");
		return false;
	}

	double inverseDeterminant = 1.0 / determinant;

	dvec3 tvec = ray.origin - a;

	double u = glm::dot(tvec, pvec) * inverseDeterminant;

	if (u < 0.0 || u > 1.0) {
		if (prd) printf("triangle no hit 2\n");
		return false;
	}

	dvec3 qvec = glm::cross(tvec, E1);
	double v = glm::dot(ray.direction, qvec) * inverseDeterminant;

	if (v < 0.0 || u + v > 1.0) {
		if (prd) printf("triangle no hit 3\n");
		return false;
	}


	hit.bary = dvec3(1.0f - (u + v), u, v);

	if(prd)printf("bary: %s\n", glm::to_string(hit.bary).c_str());

	hit.position = (a * hit.bary.x) + (b * hit.bary.y) + (c * hit.bary.z);


	hit.depth = glm::distance(hit.position, ray.origin);

	dvec3 badNormal = (this->vertices[0].normal * hit.bary.x) + (this->vertices[1].normal * hit.bary.y) + (this->vertices[2].normal * hit.bary.z);


	hit.normal = badNormal;// transformNormal(badNormal, model);
	if (prd)printf("badNormal: %s\n", glm::to_string(badNormal).c_str());
	//hit.normal = transformNormal(badNormal, dmat4(1.0));



	hit.uv = mat3x2(this->vertices[0].texCoords, this->vertices[1].texCoords, this->vertices[2].texCoords) * hit.bary;
	hit.shape = (Triangle*)this;



	//barycentric.z = glm::dot(tri.E2, qvec);

	return true;
}
void Triangle::redoAABB(double currentTime) {
	//TODO: at some point in the future this needs to take something into account

	dvec3 a = this->vertices[0].position;
	dvec3 b = this->vertices[1].position;
	dvec3 c = this->vertices[2].position;
	dvec3 minB = glm::min(a, glm::min(b, c)) - dvec3(glm::epsilon<double>());
	dvec3 maxB = glm::max(a, glm::max(b, c)) + dvec3(glm::epsilon<double>());
	this->boundingBox = AABB(minB, maxB);
}