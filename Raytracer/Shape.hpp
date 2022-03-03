#ifndef SHAPE_H
#define SHAPE_H

#include <optional>
#include <glm/glm.hpp>
#include "Ray.hpp"
#include "Material.hpp"
#include "CoordinateHelpers.hpp"

using namespace std;
using namespace glm;

struct AABB {
	dvec3 min;
	dvec3 max;
	dvec3 mid;
	AABB() {
		this->min = dvec3(0.0);
		this->max = dvec3(0.0);
	}
	AABB(dvec3 min, dvec3 max) {
		this->min = min;
		this->max = max;
		this->mid = (min + max) / 2.0;
	}
	void encompass(AABB other) {
		this->max = glm::max(this->max, other.max);
		this->min = glm::min(this->min, other.min);
	}
	bool rayAABB(const Ray& ray, optional<dvec3>& enter, optional<dvec3>& exit)const;
	bool rayAABB(const Ray& ray)const;
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

enum Side {
	RIGHT,
	LEFT,
	MIDDLE
};
class Shape
{
public:
	Shape(AABB boundingBox, const Material& material = Material("Bug"), const dmat4& model = dmat4(1.0));
	~Shape();
	bool rayAABB(const Ray& ray)const;
	virtual void redoAABB(const double currentTime) = 0;
	virtual bool rayHit(const Ray& ray, HitResult& hit, double currentTime)const = 0;
	dmat4 model;
	Material mat;
	AABB boundingBox;
protected:
private:
};











#endif