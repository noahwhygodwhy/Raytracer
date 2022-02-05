#ifndef RAY_H
#define RAY_H

#include "glm/glm.hpp"

using namespace std;
using namespace glm;

class Ray
{
public:
	dvec3 direction;
	dvec3 origin;
	Ray(dvec3 origin, dvec3 direction);
	~Ray();


private:

};

Ray::Ray(dvec3 origin, dvec3 direction)
{
	this->direction = direction;
	this->origin = origin;
}

Ray::~Ray()
{
}
#endif