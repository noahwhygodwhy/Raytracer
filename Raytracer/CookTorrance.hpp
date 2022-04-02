#ifndef COOKTORRANCE_H
#define COOKTORRANCE_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include "Ray.hpp"

#include "Material.hpp"
#include "Texture.hpp"
#include "Procedural.hpp"
#include "Shape.hpp"



using namespace std;
using namespace glm;


/*dvec3 CookTorance(
	const Ray& incomingRay,
	const Ray& outgoingRay,
	const HitResult& minRayResult,
	dvec3 downstreamRadiance,
	Material* mat,
	double prevIOR,
	double& kS);*/

dvec3 CookTorance(
	const Ray& incomingRay,
	const Ray& outgoingRay,
	const HitResult& minRayResult,
	dvec3 downstreamRadiance,
	Material* mat,
	double prevIOR,
	const dvec3& F0,
	dvec3& kS);
#endif
