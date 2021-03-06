#include "Shape.hpp"

using namespace glm;
using namespace std;



#define DIMS 3
//returns both the hit boolean, and the possibly 2 possible hit coords 
//if a enter or exit doesn't exist, it'll be NULL
//guided by https://medium.com/@bromanz/another-view-on-the-classic-ray-aabb-intersection-algorithm-for-bvh-traversal-41125138b525
bool AABB::rayAABB(const Ray& ray, optional<dvec3>& enter, optional<dvec3>& exit)const {

	dvec3 invD = 1.0 / ray.direction;
	dvec3 t0 = (this->min - ray.origin) * invD;
	dvec3 t1 = (this->max - ray.origin) * invD;

	dvec3 tSmaller = glm::min(t0, t1);
	dvec3 tBigger = glm::max(t0, t1);

	double tMin = -INFINITY;
	double tMax = INFINITY;
	tMin = glm::max(tMin, glm::max(tSmaller[0], glm::max(tSmaller[1], tSmaller[2])));
	tMax = glm::min(tMax, glm::min(tBigger[0], glm::min(tBigger[1], tBigger[2])));

	bool hit = tMin < tMax && tMax >= 0.0;

	

	enter = ray.origin + (ray.direction * tMin);
	exit = ray.origin + (ray.direction * tMax);

	return hit;
}

//pure boolean hit function
bool AABB::rayAABB(const Ray& ray)const {
	/*Side sides[DIMS]; //quadrants
	double planes[DIMS]; //candiatePlane
	double maxT[DIMS];
	bool inside = true;

	for (uint32_t i = 0; i < DIMS; i++) {
		if (ray.origin[i] < this->min[i]) {
			sides[i] = LEFT;
			planes[i] = this->min[i];
			inside = false;
		}
		else if (ray.origin[i] > this->max[i]) {
			sides[i] = RIGHT;
			planes[i] = this->max[i];
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
			if (coord[i] < this->min[i] || coord[i] > this->max[i])
				return false;
		}
		else {
			coord[i] = planes[i];
		}
	return true;*/


	//printf("\n\nray aabb\n");



	dvec3 invD = 1.0 / ray.direction;
	dvec3 t0 = (this->min - ray.origin) * invD;
	dvec3 t1 = (this->max - ray.origin) * invD;

	dvec3 tSmaller = glm::min(t0, t1);
	dvec3 tBigger = glm::max(t0, t1);

	double tMin = -INFINITY;
	double tMax = INFINITY;
	tMin = glm::max(tMin, glm::max(tSmaller[0], glm::max(tSmaller[1], tSmaller[2])));
	tMax = glm::min(tMax, glm::min(tBigger[0], glm::min(tBigger[1], tBigger[2])));

	bool hit = tMin <= tMax&& tMax >= 0.0;

	//printf("tmin: %f, tmax: %f\n", tMin, tMax);

	return hit;





}


Shape::Shape(AABB boundingBox, Material* material, const dmat4& model)
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
//based on the fast ray-box intersection chapter by andrew woo from Graphics Gems
//calculates hits between a ray and an axis aligned bounding box
bool Shape::rayAABB(const Ray& ray) const {
	return this->boundingBox.rayAABB(ray);
	//	printf("you still need to make it adjust for the transform mats");
	//	exit(0);
	
}

