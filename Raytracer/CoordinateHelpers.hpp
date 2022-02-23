#ifndef COORDINATE_HELPERS_H
#define COORDINATE_HELPERS_H

#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/glm.hpp>
using namespace std;
using namespace glm;


//dvec3 transformNormal(const dvec3& x, const dmat4& model, const dmat4& view);
dvec3 transformNormal(const dvec3& x, const dmat4& model);
dvec3 transformPos(const dvec3& x, const dmat4& model);

void solveQuadratic(double a, double b, double c, double& A1, double& A2);

#endif;

