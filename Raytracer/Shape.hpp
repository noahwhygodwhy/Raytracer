#ifndef SHAPE_H
#define SHAPE_H

#include <glm/glm.hpp>
#include "Ray.hpp"
#include "Material.hpp"
#include "CoordinateHelpers.hpp"

using namespace std;
using namespace glm;


struct AABB {
	dvec3 min;
	dvec3 max;
};

class Shape;

struct HitResult {
	dvec3 position;
	dvec3 normal;
	dvec3 tangent;
	dvec3 bitangent;
	dvec2 uv;
	dvec3 bary;//only if it's a triangle does this get used
	double depth;
	Shape* shape;
};

class Shape
{
public:
	Shape(AABB boundingBox, const Material& material, const dmat4& model);
	~Shape();
	bool rayAABB(const Ray& ray, const mat4& view);
	virtual bool rayHit(const Ray& ray, HitResult& hit, const dmat4& view, double currentTime)const = 0;
	virtual dvec3 getColor(const HitResult& hit)const = 0;
	dmat4 model;
	Material mat;
private:
	AABB boundingBox;
	
};

Shape::Shape(AABB boundingBox,const Material& material = materials.at("Bug"), const dmat4& model = dmat4(1.0))
{
	this->model = model;
	this->mat = material;
	this->boundingBox = boundingBox;
}

Shape::~Shape()
{
}

enum Side {
	RIGHT,
	LEFT,
	MIDDLE
};




//TODO: make sure this is the same as 
//https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms
#define DIMS 3
//based on the fast ray-box intersection chapter by andrew woo from Graphics Gems
//calculates hits between a ray and an axis aligned bounding box
bool Shape::rayAABB(const Ray& ray, const mat4& view) {
	printf("you still need to make it adjust for the transform mats");
	exit(0);
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




#endif