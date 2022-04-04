//#include "Triangle.hpp"
//
//
//using namespace std;
//using namespace glm;
//
//
//bool Triangle::rayHit(const Ray& ray, HitResult& hit, double currentTime) const {
//
//	
//
//
//
//	//TODO: should really be done every frame, not every frame*ray
//
//	constexpr double epsilon = glm::epsilon<double>();
//	dvec3 a = this->vertices[0].position;
//	dvec3 b = this->vertices[1].position;
//	dvec3 c = this->vertices[2].position;
//
//	dvec3 E1 = b - a;
//	dvec3 E2 = c - a;
//
//	dvec3 pvec = glm::cross(ray.direction, E2);
//	double determinant = glm::dot(E1, pvec);
//
//	if (determinant < 0.0) {//TODO: i fucked something up
//		return false;
//	}
//
//	double inverseDeterminant = 1.0 / determinant;
//
//	dvec3 tvec = ray.origin - a;
//
//	double u = glm::dot(tvec, pvec) * inverseDeterminant;
//
//	if (u < 0.0 || u > 1.0) {
//		return false;
//	}
//
//	dvec3 qvec = glm::cross(tvec, E1);
//	double v = glm::dot(ray.direction, qvec) * inverseDeterminant;
//
//	if (v < 0.0 || u + v > 1.0) {
//		return false;
//	}
//
//
//	hit.bary = dvec3(1.0f - (u + v), u, v);
//
//
//	hit.position = (a * hit.bary.x) + (b * hit.bary.y) + (c * hit.bary.z);
//
//
//	hit.depth = glm::distance(hit.position, ray.origin);
//
//	dvec3 badNormal = (this->vertices[0].normal * hit.bary.x) + (this->vertices[1].normal * hit.bary.y) + (this->vertices[2].normal * hit.bary.z);
//
//
//	hit.normal = badNormal;// transformNormal(badNormal, model);
//
//
//
//	hit.uv = mat3x2(this->vertices[0].texCoords, this->vertices[1].texCoords, this->vertices[2].texCoords) * hit.bary;
//	hit.shape = (Triangle*)this;
//
//
//
//	//barycentric.z = glm::dot(tri.E2, qvec);
//
//	return true;
//}
//void Triangle::redoAABB(double currentTime) {
//	//TODO: at some point in the future this needs to take something into account
//
//	dvec3 a = this->vertices[0].position;
//	dvec3 b = this->vertices[1].position;
//	dvec3 c = this->vertices[2].position;
//	dvec3 minB = glm::min(a, glm::min(b, c)) - dvec3(glm::epsilon<double>());
//	dvec3 maxB = glm::max(a, glm::max(b, c)) + dvec3(glm::epsilon<double>());
//	this->boundingBox = AABB(minB, maxB);
//}