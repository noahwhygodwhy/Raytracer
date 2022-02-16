#ifndef COORDINATE_HELPERS_H
#define COORDINATE_HELPERS_H

#include <glm/glm.hpp>
using namespace std;
using namespace glm;


dvec3 transformNormal(const dvec3& x, const dmat4& model) {

	dmat4 normalMat = glm::inverse(glm::transpose(model));
	return normalMat * dvec4(x, 1.0);
}
dvec3 transformPos(const dvec3& x, const dmat4& model, const dmat4& view) {
	dmat4 mvp = view * model;
	dvec4 newPos = mvp * dvec4(x, 1.0);
	return (newPos / newPos.w);
}

void solveQuadratic(double a, double b, double c, double& A1, double& A2) {
	double discriminant = (b * b) - (4.0 * a * c);
	if (discriminant < 0.0) {
		A1 = NAN;
		A2 = NAN;
		return;
	}
	A1 = (-b + glm::sqrt(discriminant)) / (2.0 * a);
	A2 = (-b - glm::sqrt(discriminant)) / (2.0 * a);
	return;
}

#endif;
