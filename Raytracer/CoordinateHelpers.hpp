#ifndef COORDINATE_HELPERS_H
#define COORDINATE_HELPERS_H



dvec3 transformNormal(const dvec3& x, const dmat4& model) {

	dmat4 normalMat = glm::inverse(glm::transpose(model));
	return normalMat * dvec4(x, 1.0);
}
dvec3 transformPos(const dvec3& x, const dmat4& model, const dmat4& view) {
	dmat4 mvp = view * model;
	dvec4 newPos = mvp * dvec4(x, 1.0f);
	return newPos / newPos.w;
}



#endif;
