#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.hpp"
#include "Animatable.hpp"


using namespace std;
using namespace glm;

class Sphere : public Shape, public Animatable
{
public:
	Sphere(const dvec3& origin, double radius, const Material& mat, dmat4(*movementFunction)(double currentTime));
	~Sphere();
	bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const;
	dvec3 getColor(const HitResult& hit)const;
	double radius;
	dvec3 origin;
	bool pointInSphere(const HitResult& hit, const dmat4& view, double currentTime);
private:
	double r2;
};

Sphere::Sphere(const dvec3& origin, double radius, const Material& mat, dmat4(*movementFunction)(double currentTime))
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

bool Sphere::pointInSphere(const HitResult& hit, const dmat4& view, double currentTime) {
	dvec3 transformedOrigin = transformPos(this->origin, this->getModel(currentTime), view);

	double bias = 1e-4;
	dvec3 biasedPoint = hit.position - (hit.normal * bias);
	
	return glm::distance(hit.position, transformedOrigin) <= this->radius;
}
dvec3 Sphere::getColor(const HitResult& hit) const {
	//return glm::abs(hit.normal);
	return vec3(1.0f, 1.0f, 1.0f);
}
//based on the geomtric sphere/ray intersection in 
//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
bool Sphere::rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime) const {

	constexpr double epsilon = glm::epsilon<double>();
	

	dmat4 currModel = this->getModel(currentTime);

	dvec3 transOrigin = transformPos(this->origin, view, currModel);

	/*dvec3 L = transOrigin - ray.origin;
	double tca = glm::dot(L, ray.direction);
	double d2 = glm::dot(L, L) - (tca * tca);
	if (d2 > this->r2) { 
		return false; 
	}
	double thc = glm::sqrt(this->r2 - d2);
	double t0 = tca - thc;
	double t1 = tca + thc;*/

	dvec3 L = ray.origin - transOrigin;
	double a = glm::dot(ray.direction, ray.direction);
	double b = 2.0 * glm::dot(ray.direction, L);
	double c = glm::dot(L, L) - this->r2;
	double t0, t1;

	double discr = (b * b) - (4.0 * a * c);
	if (discr < 0.0) { 
		return false; 
	}
	else if (discr == 0.0) { 
		t0 = t1 = -0.5 * b / a;
	}
	else {
		double q = (b > 0.0) ?
			-0.5 * (b + glm::sqrt(discr)) :
			-0.5 * (b - glm::sqrt(discr));
		t0 = q / a;
		t1 = c / q;
	}

	//double c = L.dotProduct(L) - radius2;


	if (t0 > t1) {
		swap(t0, t1);
	}
	if (t0 < 0) {
		t0 = t1;
		if (t0 < 0) { 
			return false; 
		}
	}


	//hit.normal = transformNormal(glm::normalize(hit.position - transOrigin), this->model);
	hit.position = ray.origin + (ray.direction * t0);
	//hit.normal = glm::normalize(hit.position - transOrigin);
	hit.normal = transformNormal(glm::normalize(hit.position - transOrigin), currModel);
	constexpr double spherepi = glm::pi<double>();
	hit.uv = dvec2(std::atan2(hit.normal.x, hit.normal.z) / (2 * spherepi) + 0.5, hit.normal.y * 0.5 + 0.5);

	hit.bary = dvec3(0.0);
	hit.depth = t0;
	hit.shape = (Sphere*)this;

	return true;
}
#endif;