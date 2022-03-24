#include "PointLight.hpp"

using namespace std;
using namespace glm;
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

double PointLight::getDistance(const Ray& ray) const {
	dvec3 transposition = transformPos(this->position, mat4(1));
	return glm::distance(transposition, ray.origin);
}

Ray PointLight::getRay(const dvec3& rayOrigin) const {
	//dvec3 transposition = transformPos(this->position, mat4(1));
	return Ray(rayOrigin, glm::normalize(this->position - rayOrigin));
}

double PointLight::getAttenuation(double distance)const {
	return (this->constant + (this->linear * distance) + (this->quadratic * (distance * distance)));
}