#ifndef RAY_H
#define RAY_H

#include "glm/glm.hpp"

using namespace std;
using namespace glm;

class Ray
{
public:
	vec3 direction;
	vec3 origin;
	vec3 inverseDirection;
	Ray(vec3 origin, vec3 direction);
	~Ray();


private:

};

Ray::Ray(vec3 origin, vec3 direction)
{
	this->direction = direction;
	this->origin = origin;
	this->inverseDirection = -direction;
}

Ray::~Ray()
{
}
#endif