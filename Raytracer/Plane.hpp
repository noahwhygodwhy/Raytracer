#ifndef PLANE_H
#define PLANE_H

#include <glm/glm.hpp>

#include "Shape.hpp"

using namespace std;
using namespace glm;

class Plane : public Shape
{
public:
	Plane(const dvec3& a, const dvec3& b);
	~Plane();
	bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view)const;
	dvec3 getColor(const HitResult& hit)const;
private:
	dvec3 a, b;
};

Plane::Plane(const dvec3& a, const dvec3& b) : Shape(AABB(a, b))
{
	this->a = a;
	this->b = b;
}

Plane::~Plane()
{
}

bool Plane::rayHit(const Ray& ray, HitResult& hit, const dmat4& view) const {
	return false;
}


dvec3 Plane::getColor(const HitResult& hit) const {
	return dvec3(1.0);
}




#endif