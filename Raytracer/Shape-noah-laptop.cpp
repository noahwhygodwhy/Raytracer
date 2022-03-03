#include "Shape.hpp"

using namespace glm;
using namespace std;



void AABB::rayAABB(const Ray& ray, dvec3& enter, dvec3& exit) {

}

Shape::Shape(AABB boundingBox, const Material& material, const dmat4& model)
{
	this->model = model;
	this->mat = material;
	this->boundingBox = boundingBox;
}

Shape::~Shape()
{
}


//TODO: make sure this is the same as 
//https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
#define DIMS 3
//based on the fast ray-box intersection chapter by andrew woo from Graphics Gems
//calculates hits between a ray and an axis aligned bounding box
bool Shape::rayAABB(const Ray& ray) const {
	//	printf("you still need to make it adjust for the transform mats");
	//	exit(0);
	Side sides[DIMS]; //quadrants
	double planes[DIMS]; //candiatePlane
	double maxT[DIMS];
	bool inside = true;

	for (uint32_t i = 0; i < DIMS; i++) {
		if (ray.origin[i] < this->boundingBox.min[i]) {
			sides[i] = LEFT;
			planes[i] = this->boundingBox.min[i];
			inside = false;
		}
		else if (ray.origin[i] > this->boundingBox.max[i]) {
			sides[i] = RIGHT;
			planes[i] = this->boundingBox.max[i];
			inside = false;
		}
		else {
			sides[i] = MIDDLE;
		}
	}
	if (inside) {
		return true;
	}
	for (uint32_t i = 0; i < DIMS; i++)
	{
		if (sides[i] != MIDDLE && ray.direction[i] != 0.0) {
			maxT[i] = (planes[i] - ray.origin[i]) / ray.direction[i];
		}
		else {
			maxT[i] = -1.0;
		}
	}

	uint32_t intersectionDim = 0;
	for (uint32_t i = 1; i < DIMS; i++) {
		if (maxT[intersectionDim] < maxT[i])
		{
			intersectionDim = i;
		}
	}

	if (maxT[intersectionDim] < 0.0f) {
		return (false);
	}
	dvec3 coord;

	for (uint32_t i = 0; i < DIMS; i++)
		if (intersectionDim != i) {
			coord[i] = ray.origin[i] + maxT[intersectionDim] * ray.direction[i];
			if (coord[i] < this->boundingBox.min[i] || coord[i] > this->boundingBox.max[i])
				return false;
		}
		else {
			coord[i] = planes[i];
		}
	return true;
}

