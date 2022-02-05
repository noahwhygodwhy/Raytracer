#ifndef DIRECTIONAL_LIGHT_H
#define DIRECTIONAL_LIGHT_H

#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class DirectionalLight
{
public:
	DirectionalLight(vec3 direction, vec3 color);
	~DirectionalLight();

private:

};

DirectionalLight::DirectionalLight()
{
}

DirectionalLight::~DirectionalLight()
{
}
#endif;