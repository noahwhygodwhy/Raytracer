#ifndef DIRECTIONAL_LIGHT_H
#define DIRECTIONAL_LIGHT_H

#include <glm/glm.hpp>
#include "Ray.hpp"
using namespace std;
using namespace glm;



constexpr double bigNumberButNotInfinity = 10000000000000.0;

class DirectionalLight : public Light
{
public:
	DirectionalLight(dvec3 direction, dvec3 color);
	~DirectionalLight();
	double getDistance(const dvec3& rayOrigin) const;
	Ray getRay(const dvec3& rayOrigin) const;
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

double DirectionalLight::getDistance(const dvec3& rayOrigin) const {
	return bigNumberButNotInfinity;
}

Ray DirectionalLight::getRay(const dvec3& rayOrigin) const {


	//dvec3 alteredRayOrigin = transformPos(rayOrigin, mat4(1.0), view);

	dvec3 alteredReverseDirection = dvec4(this->reverseDirection, 1.0f);

	dvec3 transposition = transformPos(this->reverseDirection*bigNumberButNotInfinity, mat4(1));

	return Ray(rayOrigin, glm::normalize(transposition - rayOrigin));
}
double DirectionalLight::getAttenuation(double distance)const {
	return 1.0;
}

#endif;