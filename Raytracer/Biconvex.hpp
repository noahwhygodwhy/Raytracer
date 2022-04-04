//#ifndef BICONVEX_H
//#define BICONVEX_H
//
//#include <glm/glm.hpp>
//#include "Shape.hpp"
//#include "Animatable.hpp"
//#include "Sphere.hpp"
//
//using namespace std;
//using namespace glm;
//
//
//AABB makeBoundingBox(double sphereRadius, double focalLength, double ior);
//double solveLensForD(double r, double f, double n);
//
//class Biconvex : public Shape, public Animatable
//{
//public:
//
//	Biconvex(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength, dmat4(*movementFunction)(double currentTime))
//		: Shape(this->setupSpheres(origin, forward, lensRadius, focalLength))
//		, Animatable(movementFunction)
//	{
//		//nothing here cause it all really happens in setupSpheres
//	}	~Biconvex();
//	bool rayHit(const Ray& ray, HitResult& hit, double currentTime)const;
//private:
//	Sphere* s1;
//	Sphere* s2;
//	double d;
//	double focal;
//	dvec3 origin;
//	dvec3 forward;
//	AABB setupSpheres(const dvec3 origin, dvec3 forward, double lensRadius, double focalLength);
//};
//
//
//
//
//
//#endif