#include "Procedural.hpp"


dvec3 ryCheckers10x10(dvec2 uv) {
	dvec2 squares(1000.0, 1000.0);
	ivec2 flatUV = glm::floor(uv * squares);
	if ((flatUV.x + flatUV.y) % 2 == 0) {
		return dvec3(1.0, 0.0, 0.0);

	}
	else {
		return dvec3(1.0, 1.0, 0.0);
	}
}

double doubleCheckers10x10(dvec2 uv) {
	dvec2 squares(10.0, 10.0);
	ivec2 flatUV = glm::floor(uv * squares);
	if ((flatUV.x + flatUV.y) % 2 == 0) {
		return 1.0;
	}
	else {
		return 0.0;
	}
}