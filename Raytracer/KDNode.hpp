#ifndef KDNODE_H
#define KDNODE_H

#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <execution>

#include "Triangle.hpp"


#define OBJECTS_PER_LEAF 5

enum Axis {
	X,
	Y,
	Z
};

using namespace std;
using namespace glm;




class KDNode
{
public:
	KDNode();
	~KDNode();

private:

};

KDNode::KDNode()
{
}

KDNode::~KDNode()
{
}



class KDBranch:KDNode
{
public:
	KDBranch(vector<Triangle> triangles, Axis axis);
	~KDBranch();
	//vector<Triangle> intersects;
	KDNode* lesser;
	KDNode* greater;
	Axis axis;
	double value;

private:

};


const char* toBool(bool x) {
	return x ? "True" : "False";
}




//TODO: perhaps consider instead of having an intersects, just give those that
//intersect to both left and right, to be sorted at a future time, with only leaf nodes
//having any triangles at all.

KDBranch::KDBranch(vector<Triangle> triangles, Axis axis = X)
{

	this->axis = axis;
	auto compTri = [axis](const Triangle& a, const Triangle& b)
	{
		return int32_t(a.middleBounding[axis] < b.middleBounding[axis]?1:0);
	};
	std::sort(execution::par_unseq, triangles.begin(), triangles.end(), compTri);
	int32_t midIndex = triangles.size() / 2;
	
	this->value = triangles.at(midIndex).middleBounding[axis];
	int32_t upperIndex = midIndex;
	int32_t lowerIndex = midIndex;
	while (upperIndex < triangles.size()-1 && triangles.at(upperIndex).minBounding[axis] < this->value) {
		upperIndex++;
	}
	while (lowerIndex > 0 && triangles.at(lowerIndex).maxBounding[axis] > this->value) {
		lowerIndex--;
	}

	Axis nextAxis = Axis((axis + 1) % 3);
	//this->intersects = vector<Triangle>(triangles.begin() + lowerIndex, triangles.begin() + upperIndex+1);

	if (lowerIndex < 1) {
		this->lesser = NULL;
	} else {
		this->lesser = new KDBranch(vector<Triangle>(triangles.begin() + 0, triangles.begin() + lowerIndex), nextAxis);
	}
	if (upperIndex > (int32_t(triangles.size())-2)) {
		this->greater = NULL;
	} else {
		this->greater = new KDBranch(vector<Triangle>(triangles.begin() + upperIndex+1, triangles.end()), nextAxis);
	}
}

KDBranch::~KDBranch()
{
}


enum KDCase {
	POSITIVE,
	NEGATIVE,
	POSTONEG,
	NEGTOPOS
};



KDCase rayHit(const KDBranch* kdn, const Ray& ray) {
	double rayvalue = ray.origin[kdn->axis];
	double raydir = ray.direction[kdn->axis];
	double epsilon = glm::epsilon<double>();
	if (rayvalue > kdn->value) {
		if (raydir <= epsilon) {
			return POSTONEG;
		}
		return POSITIVE;
	}
	else {
		if (raydir >= -epsilon) {
			return NEGTOPOS;
		}
		return NEGATIVE;
	}
}





Triangle* evaluateIntersects(const KDBranch* kdn, const Ray& ray, vec3& barycentric, double& interDepth) {
	
	vec3 hitPoint;
	Triangle* closest = NULL;
	double closestDepth = INFINITY;
	
	for (const Triangle& tri : kdn->intersects) {
		if (rayHit(ray, tri.minBounding, tri.maxBounding)) {
			if (rayHit(ray, tri, barycentric, hitPoint)) {
				double depth = glm::length2(hitPoint-ray.origin);
				if (depth < closestDepth) {
					closest = (Triangle*) & tri;
					closestDepth = depth;
				}
			}
		}
	}
	interDepth = closestDepth;
	return closest;
}


Triangle* rayHit(const KDBranch* kdn, const Ray& ray, vec3& barycentric, double& depthOut) {
	




	double trueDepth = INFINITY;
	Triangle* trueT = NULL;

	Triangle* posT = NULL;
	double posDepth = INFINITY;

	Triangle* interT = NULL;
	double interDepth = INFINITY;

	Triangle* negT = NULL;
	double negDepth = INFINITY;

	if (kdn->greater != NULL) {
		posT = rayHit(kdn->greater, ray, barycentric, posDepth);
	}
	interT = evaluateIntersects(kdn, ray, barycentric, interDepth);
	if (kdn->lesser != NULL) {
		negT = rayHit(kdn->lesser, ray, barycentric, negDepth);
	}

	if (posDepth < trueDepth && posT != NULL) {
		trueDepth = posDepth;
		trueT = posT;
	}
	if (negDepth < trueDepth && negT != NULL) {
		trueDepth = negDepth;
		trueT = negT;
	}
	if (interDepth < trueDepth && interT != NULL) {
		trueDepth = interDepth;
		trueT = interT;
	}
	depthOut = trueDepth;
	return trueT;


	/*
	switch (rayHit(kdn, ray)) {
		case POSITIVE:
			if (kdn->greater != NULL) {
				t = rayHit(kdn->greater, ray, barycentric);
			}
			if (t != NULL) {
				t = evaluateIntersects(kdn, ray, barycentric);
			}
			//rayHit() but for the positive side of kdn,
			//then evaluate the intersect ones for kdn
			break;
		case NEGATIVE:
			if (kdn->lesser != NULL) {
				t = rayHit(kdn->lesser, ray, barycentric);
			}
			if (t == NULL) {
				t = evaluateIntersects(kdn, ray, barycentric);
			}
			//rayHit() but for the negative side of kdn,
			//then evaluate the intersect ones for kdn
			break;
		case POSTONEG:
			if (kdn->greater != NULL) {
				t = rayHit(kdn->greater, ray, barycentric);
			}
			if (t == NULL) {
				t = evaluateIntersects(kdn, ray, barycentric);
			}
			if (t == NULL && kdn->lesser != NULL) {
				t = rayHit(kdn->lesser, ray, barycentric);
			}
			//rayHit() but for the positive side of kdn,
			//then evaluate the intersect ones for kdn
			//rayHit() but for the negative side of kdn,
			break;
		case NEGTOPOS:
			if (kdn->lesser != NULL) {
				t = rayHit(kdn->lesser, ray, barycentric);
			}
			if (t == NULL) {
				t = evaluateIntersects(kdn, ray, barycentric);
			}
			if (t == NULL && kdn->greater != NULL) {
				t = rayHit(kdn->greater, ray, barycentric);
			}
			//rayHit() but for the negative side of kdn,
			//then evaluate the intersect ones for kdn
			//rayHit() but for the positive side of kdn,
			break;
		default:
			printf("\n\nBad KDCase\n\n");
	}
	return t;*/
}

#endif;
