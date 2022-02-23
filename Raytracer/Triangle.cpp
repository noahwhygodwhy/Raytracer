#include "Triangle.hpp"


using namespace std;
using namespace glm;


bool Triangle::rayHit(const Ray& ray, HitResult& hit, double currentTime) const {


	//TODO: should really be done every frame, not every frame*ray
	constexpr double epsilon = glm::epsilon<double>();
	dvec3 a = transformPos(this->vertices[0].position, dmat4(1.0));
	dvec3 b = transformPos(this->vertices[1].position, dmat4(1.0));
	dvec3 c = transformPos(this->vertices[2].position, dmat4(1.0));

	dvec3 E1 = b - a;
	dvec3 E2 = c - a;

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

	if (prd)printf("hit at %s\n", glm::to_string(hit.position).c_str());

	hit.depth = glm::length(hit.position - ray.origin);
	//hit.normal = glm::normalize(transformedNormal);




	dvec3 badNormal = (this->vertices[0].normal * hit.bary.x) + (this->vertices[1].normal * hit.bary.y) + (this->vertices[2].normal * hit.bary.z);

	if (prd)printf("vert 1 normal: %s\n", glm::to_string(this->vertices[0].normal).c_str());
	if (prd)printf("vert 2 normal: %s\n", glm::to_string(this->vertices[1].normal).c_str());
	if (prd)printf("vert 3 normal: %s\n", glm::to_string(this->vertices[2].normal).c_str());
	if (prd)printf("bary: %s\n", glm::to_string(hit.bary).c_str());


	hit.normal = transformNormal(badNormal, model);
	if (prd)printf("normal: %s\n", glm::to_string(hit.normal).c_str());

	//hit.normal = transformNormal(badNormal, dmat4(1.0));



	hit.uv = mat3x2(this->vertices[0].texCoords, this->vertices[1].texCoords, this->vertices[2].texCoords) * hit.bary;
	hit.shape = (Triangle*)this;



	//barycentric.z = glm::dot(tri.E2, qvec);

	return true;
}
void Triangle::redoAABB(double currentTime) {

	//TODO: at some point in the future this needs to take something into account

	dvec3 a = transformPos(this->vertices[0].position, dmat4(1.0));
	dvec3 b = transformPos(this->vertices[1].position, dmat4(1.0));
	dvec3 c = transformPos(this->vertices[2].position, dmat4(1.0));

	this->boundingBox = AABB(glm::min(a, glm::min(b, c)), glm::max(a, glm::max(b, c)));
}