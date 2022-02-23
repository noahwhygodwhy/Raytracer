#include "Sphere.hpp"

using namespace std;
using namespace glm;

Sphere::Sphere(const dvec3& origin, double radius, Material mat, dmat4(*movementFunction)(double currentTime))
	: Shape(AABB(origin - dvec3(radius), origin + dvec3(radius)), mat)
	, Animatable(movementFunction)
{
	this->origin = origin;
	this->radius = radius;
	this->r2 = radius * radius;

}

Sphere::~Sphere()
{
}

bool Sphere::pointInSphere(const HitResult& hit, double currentTime) {
	dvec3 transformedOrigin = transformPos(this->origin, this->getModel(currentTime));

	double bias = 1e-4;
	dvec3 biasedPoint = hit.position - (hit.normal * bias);

	return glm::distance(hit.position, transformedOrigin) <= this->radius;
}




void Sphere::redoAABB(double currentTime) {
	dvec3 transOrigin = transformPos(this->origin, this->getModel(currentTime));
	this->boundingBox = AABB(transOrigin - dvec3(radius), transOrigin + dvec3(radius));
}



//based on the geomtric sphere/ray intersection in 
//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
bool Sphere::rayHit(const Ray& ray, HitResult& hit, double currentTime) const {

	
	constexpr double epsilon = glm::epsilon<double>();


	dmat4 currModel = this->getModel(currentTime);
	dvec3 transOrigin = transformPos(this->origin, currModel);


	dvec3 L = ray.origin - transOrigin;
	double a = glm::dot(ray.direction, ray.direction);
	double b = 2.0 * glm::dot(ray.direction, L);
	double c = glm::dot(L, L) - this->r2;
	double t0, t1;
	solveQuadratic(a, b, c, t0, t1);

	if (t0 == NAN && t1 == NAN) {
		return false;
	}

	if (t0 > t1) {
		swap(t0, t1);
	}
	if (t0 < 0) {
		t0 = t1;
		if (t0 < 0) {
			return false;
		}
	}


	hit.position = ray.origin + (ray.direction * t0);
	hit.normal = transformNormal(glm::normalize(hit.position - transOrigin), currModel);
	hit.uv = dvec2(std::atan2(hit.normal.x, hit.normal.z) / (2 * glm::pi<double>()) + 0.5, hit.normal.y * 0.5 + 0.5);
	hit.bary = dvec3(0.0);
	hit.depth = t0;
	hit.shape = (Sphere*)this;

	return true;
}