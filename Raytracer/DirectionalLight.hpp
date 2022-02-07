#ifndef DIRECTIONAL_LIGHT_H
#define DIRECTIONAL_LIGHT_H

#include <glm/glm.hpp>
#include "Ray.hpp"
using namespace std;
using namespace glm;

class DirectionalLight : public Light
{
public:
	DirectionalLight(dvec3 direction, dvec3 color);
	~DirectionalLight();
	double getDistance(const dvec3& rayOrigin, const dmat4& view) const;
	Ray getRay(const dvec3& rayOrigin, const dmat4& view) const;
	double getAttenuation(double distance)const;
private:
	dvec3 direction;
	dvec3 reverseDirection;
};

DirectionalLight::DirectionalLight(dvec3 direction, dvec3 color) : Light(color)
{
	this->direction = glm::normalize(direction);
	this->reverseDirection = -this->direction;
}

DirectionalLight::~DirectionalLight()
{
}

double DirectionalLight::getDistance(const dvec3& rayOrigin, const dmat4& view) const {
	return 10000000000000.0;//cheating, but hey
}

Ray DirectionalLight::getRay(const dvec3& rayOrigin, const dmat4& view) const {
	return Ray(rayOrigin, this->reverseDirection);
}
double DirectionalLight::getAttenuation(double distance)const {
	return 1.0;
}

#endif;