#ifndef LIGHT_H
#define LIGHT_H

#include <glm/glm.hpp>


#include "Ray.hpp"

using namespace std;
using namespace glm;


#define GLOBAL_AMBIENT 0.1


class Light
{
public:
	Light(const dvec3& color);
	~Light();
	//virtual vec3 rayHit(const vec3& rayOrigin, const mat4& view)const = 0;
	virtual double getDistance(const dvec3& origin) const = 0;
	virtual Ray getRay(const dvec3& rayOrigin) const = 0;
	virtual double getAttenuation(double distance)const = 0;
	dvec3 color;
private:

};

Light::Light(const dvec3& color)
{
	this->color = color;
}

Light::~Light()
{
}

#endif