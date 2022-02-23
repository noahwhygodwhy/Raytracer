#ifndef SPHERE_H
#define SPHERE_H
#include "Shape.hpp"
#include "Animatable.hpp"
#include "Material.hpp"
#include "CoordinateHelpers.hpp"
#include <glm/gtc/epsilon.hpp>
#include <glm/glm.hpp>


extern bool prd;

using namespace std;
using namespace glm;

class Sphere : public Shape, public Animatable
{
public:
	Sphere(const dvec3& origin, double radius, Material, dmat4(*movementFunction)(double currentTime));
	~Sphere();
	bool rayHit(const Ray& ray, HitResult& hit, double currentTime)const;
	double radius;
	dvec3 origin;
	bool pointInSphere(const HitResult& hit, double currentTime);
	void redoAABB(double currentTime);
private:
	double r2;
};

#endif;