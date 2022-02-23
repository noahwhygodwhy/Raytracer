#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H

#include <glm/glm.hpp>

#include "Light.hpp"

using namespace std;
using namespace glm;

class PointLight : public Light
{
public:
	PointLight(const dvec3& position, const dvec3& color, const dvec3& attenuationVals);
	~PointLight();
	double getDistance(const dvec3& rayOrigin) const;
	Ray getRay(const dvec3& rayOrigin) const;
	double getAttenuation(double distance)const;
private:
	dvec3 position; 
	double constant;
	double linear;
	double quadratic;

};


PointLight::PointLight(const dvec3& position, const dvec3& color, const dvec3& attenuationVals) : Light(color)
{
	this->position = position;
	this->constant = attenuationVals.x;
	this->linear = attenuationVals.y;
	this->quadratic = attenuationVals.z;

}

PointLight::~PointLight()
{
}

/*vec3 PointLight::rayHit(vec3 rayOrigin, const mat4& view)const {

}*/
double PointLight::getDistance(const dvec3& rayOrigin) const {
	dvec3 transposition = transformPos(this->position, mat4(1));
	return glm::distance(transposition, rayOrigin);
}

Ray PointLight::getRay(const dvec3& rayOrigin) const {
	dvec3 transposition = transformPos(this->position, mat4(1));

	return Ray(rayOrigin, glm::normalize(transposition-rayOrigin));
}

double PointLight::getAttenuation(double distance)const {

	//return (((this->constant / distance) + this->linear + (this->quadratic * distance)) * distance
	//return (this->constant + ((this->linear + (this->quadratic * distance)) * distance));
	return (this->constant + (this->linear * distance) + (this->quadratic * (distance * distance)));
}

#endif;