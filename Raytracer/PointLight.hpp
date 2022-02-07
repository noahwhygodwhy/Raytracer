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
	double getDistance(const dvec3& rayOrigin, const dmat4& view) const;
	Ray getRay(const dvec3& rayOrigin, const dmat4& view) const;
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
double PointLight::getDistance(const dvec3& rayOrigin, const dmat4& view) const {
	dvec3 transposition = transformPos(this->position, mat4(1), view);
	/*printf("get distance 2 of light\nposition:%s\ntransposition: %s\nray origin: %s\ndistance2:%f\n",
		glm::to_string(this->position).c_str(),
		glm::to_string(transposition).c_str(),
		glm::to_string(rayOrigin).c_str(),
		glm::distance2(transposition, rayOrigin));*/
	//return glm::length2(transposition - rayOrigin);

	return glm::distance(transposition, rayOrigin);
}

Ray PointLight::getRay(const dvec3& rayOrigin, const dmat4& view) const {
	dvec3 transposition = transformPos(this->position, mat4(1), view);

	return Ray(rayOrigin, glm::normalize(transposition -rayOrigin));
}

double PointLight::getAttenuation(double distance)const {

	//return (((this->constant / distance) + this->linear + (this->quadratic * distance)) * distance
	return (this->constant + ((this->linear + (this->quadratic * distance)) * distance));
	return (this->constant + (this->linear * distance) + (this->quadratic * (distance * distance)));
}

#endif;