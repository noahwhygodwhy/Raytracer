#ifndef RAY_H
#define RAY_H

#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Ray
{
public:
	dvec3 direction;
	dvec3 inverseDirection;
	dvec3 origin;
	Ray(const dvec3& origin, const dvec3& direction);
	~Ray();


private:

};

#endif