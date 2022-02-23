#include "Biconvex.hpp"

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



/*Biconvex::Biconvex(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength, dmat4(*movementFunction)(double currentTime))
	: Shape(this->setupSpheres(origin, forward, lensRadius, focalLength))
	, Animatable(movementFunction)
{
	//nothing here cause it all really happens in setupSpheres
}*/

AABB Biconvex::setupSpheres(const dvec3 origin, dvec3 forward, double sphereRadius, double focalLength) {
	Material glass("Glass");
	this->d = solveLensForD(sphereRadius, focalLength, 1.54);//TODO: this can't be hardcoded, screw the focal length


	this->origin = origin;
	this->forward = glm::normalize(forward);



	dvec3 s1Origin = origin - (this->forward * (sphereRadius - d));
	dvec3 s2Origin = origin + (this->forward * (sphereRadius - d));


	this->s1 = new Sphere(s1Origin, sphereRadius, glass, noMovement);
	this->s2 = new Sphere(s2Origin, sphereRadius, glass, noMovement);

	double lbsr = sqrt((sphereRadius - (d / 2)) * (8.0 * d)) / 2.0; //lense bounding sphere radius

	//TODO: need to realculate this on rotation or view change ughhh how do i do this easily
	//need a getAABB() functio
	return AABB(origin - dvec3(lbsr), origin + dvec3(lbsr));
}



bool Biconvex::rayHit(const Ray& ray, HitResult& hit, double currentTime)const {
	HitResult s1HitRes, s2HitRes;
	bool s1Hit = this->s1->rayHit(ray, s1HitRes, currentTime);
	bool s2Hit = this->s2->rayHit(ray, s2HitRes, currentTime);



	dvec3 transformeds1Origin = transformPos(s1->origin, s1->getModel(currentTime));
	dvec3 transformeds2Origin = transformPos(s2->origin, s2->getModel(currentTime));

	//printf("\n\n============================\n");
	//printf("s1HitPos: %s\n", glm::to_string(s1HitRes.position).c_str());
	//printf("s1origin: s1Radius: %s, %d\n", glm::to_string(transformeds1Origin).c_str(), s1->radius);
	//printf("s2HitPos: %s\n", glm::to_string(s2HitRes.position).c_str());
	//printf("s2origin: s2Radius: %s, %d\n", glm::to_string(transformeds2Origin).c_str(), s2->radius);

	s1Hit = s1Hit && s2->pointInSphere(s1HitRes, currentTime);
	s2Hit = s2Hit && s1->pointInSphere(s2HitRes, currentTime);

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



Biconvex::~Biconvex()
{
}