#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.hpp"
#include "Animatable.hpp"


using namespace std;
using namespace glm;

class Sphere : public Shape, public Animatable
{
public:
	Sphere(const dvec3& origin, double radius, dmat4(*movementFunction)(double currentTime));
	~Sphere();
	bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const;
	dvec3 getColor(const HitResult& hit)const;
private:
	double radius;
	double r2;
	dvec3 origin;
};

Sphere::Sphere(const dvec3& origin, double radius, dmat4(*movementFunction)(double currentTime))
	: Shape(AABB(origin - dvec3(radius), origin + dvec3(radius)))
	, Animatable(movementFunction)
{
	this->origin = origin;
	this->radius = radius;
	this->r2 = radius * radius;
}

Sphere::~Sphere()
{
}


dvec3 Sphere::getColor(const HitResult& hit) const {
	//return glm::abs(hit.normal);
	return vec3(1.0f, 0.0f, 0.2f);
}
//based on the geomtric sphere/ray intersection in 
//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
bool Sphere::rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime) const {

	constexpr double epsilon = glm::epsilon<double>();
	
	dvec3 transOrigin = transformPos(this->origin, view, this->getModel(currentTime));

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


	hit.position = ray.origin + (ray.direction * dvec3(t0));
	//hit.normal = transformNormal(glm::normalize(hit.position - transOrigin), this->model);
	hit.normal = glm::normalize(hit.position - transOrigin);
	hit.bary = dvec3(0.0);
	hit.depth = t0;
	hit.shape = (Sphere*)this;

	return true;
}
#endif;