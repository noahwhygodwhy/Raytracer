#ifndef BICONVEX_H
#define BICONVEX_H

#include "glm/glm.hpp"
#include "Shape.hpp"
#include "Animatable.hpp"
#include "Sphere.hpp"

using namespace std;
using namespace glm;




double solveLensForD(double r, double f, double n) {

	double numerator = n * r * r;
	double denom = f * (n - 1) * (n - 1);
	return numerator / denom;

	//return (((1 / f) / (n - 1)) * (n * r * r)) / (n - 1);
}



AABB makeBoundingBox(double sphereRadius, double focalLength, double ior) {
	double d = solveLensForD(sphereRadius, focalLength, ior);
	//TODO:
	return AABB();
}

class Biconvex : public Shape, public Animatable
{
public:
	Biconvex(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength, dmat4(*movementFunction)(double currentTime));
	~Biconvex();
	bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const;
	dvec3 getColor(const HitResult& hit)const;
private:
	Sphere* s1;
	Sphere* s2;
	double d;
	double focal;
	dvec3 origin;
	dvec3 forward;
	AABB setupSpheres(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength);
};

Biconvex::Biconvex(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength, dmat4(*movementFunction)(double currentTime))
	: Shape(this->setupSpheres(origin, forward, lensRadius, focalLength))
	, Animatable(movementFunction)
{
	//nothing here cause it all really happens in setupSpheres
}

AABB Biconvex::setupSpheres(const dvec3 origin, dvec3 forward, double sphereRadius, double focalLength) {
	Material glass = materials.at("Glass");
	this->d = solveLensForD(sphereRadius, focalLength, glass.ni);

	printf("d: %f\n", d);

	this->origin = origin;
	this->forward = glm::normalize(forward);


	printf("forward: %s\n", glm::to_string(this->forward).c_str());

	dvec3 s1Origin = origin - (this->forward * (sphereRadius - d));
	dvec3 s2Origin = origin + (this->forward * (sphereRadius - d));

	printf("making sphere with origin %s\n", glm::to_string(s1Origin).c_str());
	printf("making sphere with origin %s\n", glm::to_string(s2Origin).c_str());

	this->s1 = new Sphere(s1Origin, sphereRadius, glass, noMovement);
	this->s2 = new Sphere(s2Origin, sphereRadius, glass, noMovement);

	double lbsr = sqrt((sphereRadius - (d / 2)) * (8.0 * d)) / 2.0; //lense bounding sphere radius

	//TODO: need to realculate this on rotation or view change ughhh how do i do this easily
	//need a getAABB() functio
	return AABB(origin - dvec3(lbsr), origin + dvec3(lbsr));
}



bool Biconvex::rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const {
	HitResult s1HitRes, s2HitRes;
	bool s1Hit = this->s1->rayHit(ray, s1HitRes, view, currentTime);
	bool s2Hit = this->s2->rayHit(ray, s2HitRes, view, currentTime);



	dvec3 transformeds1Origin = transformPos(s1->origin, s1->getModel(currentTime), view);
	dvec3 transformeds2Origin = transformPos(s2->origin, s2->getModel(currentTime), view);

	//printf("\n\n============================\n");
	//printf("s1HitPos: %s\n", glm::to_string(s1HitRes.position).c_str());
	//printf("s1origin: s1Radius: %s, %d\n", glm::to_string(transformeds1Origin).c_str(), s1->radius);
	//printf("s2HitPos: %s\n", glm::to_string(s2HitRes.position).c_str());
	//printf("s2origin: s2Radius: %s, %d\n", glm::to_string(transformeds2Origin).c_str(), s2->radius);

	s1Hit = s1Hit && s2->pointInSphere(s1HitRes, view, currentTime);
	s2Hit = s2Hit && s1->pointInSphere(s2HitRes, view, currentTime);

	double minDistance = INFINITY;

	if (s1Hit && s1HitRes.depth < minDistance) {
		minDistance = s1HitRes.depth;
		hit = s1HitRes;
	}
	if (s2Hit && s2HitRes.depth < minDistance) {
		hit = s2HitRes;
	}

	return s1Hit || s2Hit;
}
dvec3 Biconvex::getColor(const HitResult& hit)const {
	return dvec3(1.0);
}


Biconvex::~Biconvex()
{
}

#endif