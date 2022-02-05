//#ifndef KDNODE_H
//#define KDNODE_H
//
//#include <vector>
//#include <algorithm>
//#include <execution>
//#include <stack>
//
//
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
//#include <glm/gtx/string_cast.hpp>
//#include <glm/gtx/norm.hpp>
//
//#include "Triangle.hpp"
//
//
//
//using namespace std;
//using namespace glm;
//
//uint32_t OBJECTS_PER_LEAF = 10;
//
//enum Axis {
//	X,
//	Y,
//	Z
//};
//
//
//
//
///*class ImpKDNode
//{
//public:
//	bool traversed;
//	//uvec3 U;
//	vector<Triangle> triangles;
//	ImpKDNode() {}
//	~ImpKDNode() {}
//
//
//private:
//
//};*/
//
//
//
//
//
//
//
//class KDNode
//{
//public:
//	bool isLeaf;
//	//bool traversed;
//	//KDNode* parent;
//	//uvec3 U;
//	KDNode() {
//		//this->parent = parent;
//	}
//	~KDNode() {}
//	
//
//private:
//
//};
//
//
//
//
//class KDBranch : public KDNode
//{
//public:
//	KDBranch(vector<Triangle>* triangles, Axis axis, int layer);
//	KDBranch(KDNode* lesser, KDNode* greater, Axis axis, double value);
//	~KDBranch();
//	//vector<Triangle> intersects;
//	KDNode* lesser;
//	KDNode* greater;
//	Axis axis;
//	double value;
//
//private:
//
//};
//
//class KDLeaf : public KDNode
//{
//public:
//	KDLeaf(vector<Triangle>* tris, Axis axis):KDNode() {
//		this->tris = tris;
//		this->isLeaf = true;
//	}
//	~KDLeaf() {}
//	vector<Triangle>* tris;
//private:
//
//};
//
//
//const char* toBool(bool x) {
//	return x ? "True" : "False";
//}
//
//
//
//KDBranch::KDBranch(KDNode* lesser, KDNode* greater, Axis axis, double value) {
//	this->isLeaf = false;
//}
//
//
//KDNode* makeKDTree(vector<Triangle>* triangles, Axis axis = X, uint32_t lowest = UINT32_MAX, uint32_t count = 0, uint32_t layer = 0) {
//	auto compTri = [axis](const Triangle& a, const Triangle& b)
//	{
//		return int32_t(a.middleBounding[axis] < b.middleBounding[axis] ? 1 : 0);
//	};
//	printf("%u======================\ntriangles: %li\n", layer, triangles->size());
//	printf("lowest: %i, count: %i\n", lowest, count);
//	if (triangles->size() >= lowest) {
//		count++;
//	}
//	else {
//		count = 0;
//	}
//	if (triangles->size() < OBJECTS_PER_LEAF || count > 2) {
//		printf("making leaf\n");
//		return new KDLeaf(triangles, axis);
//	}
//	lowest = glm::min(lowest, uint32_t(triangles->size()));
//
//	std::sort(execution::par_unseq, triangles->begin(), triangles->end(), compTri);
//	int32_t midIndex = triangles->size() / 2;
//	double value = triangles->at(midIndex).middleBounding[axis];
//	int32_t upperIndex = midIndex;
//	int32_t lowerIndex = midIndex;
//
//	if (triangles->size() >= lowest) {
//		count++;
//	} else {
//		count = 0;
//	}
//	lowest = glm::min(lowest, uint32_t(triangles->size()));
//	if (count > 2) {
//
//		count = 0;
//	}
//	while (upperIndex < triangles->size() - 1 && triangles->at(upperIndex).minBounding[axis] < value) {
//		upperIndex++;
//	}
//	while (lowerIndex > 0 && triangles->at(lowerIndex).maxBounding[axis] > value) {
//		lowerIndex--;
//	}
//	ptrdiff_t x = triangles->end() - triangles->begin();
//	printf("upper: %i, lower: %i\n", upperIndex, lowerIndex);
//	printf("lower is %i to %li\n", 0, triangles->end() - (triangles->begin() + lowerIndex));
//	printf("greater is %i to %i\n", triangles->end() - (triangles->begin() + lowerIndex), triangles->end() - triangles->begin());
//	Axis nextAxis = Axis((axis + 1) % 3);
//	return new KDBranch(
//		makeKDTree(new vector<Triangle>(triangles->begin(), triangles->end() - lowerIndex), nextAxis, lowest, count, layer+1),
//		makeKDTree(new vector<Triangle>(triangles->end()-upperIndex, triangles->end()), nextAxis, lowest, count, layer+1),
//		axis,
//		value
//	);
//
//
//
//	//->lesser = new KDBranch(new vector<Triangle>(triangles->begin() + 0, triangles->begin() + lowerIndex), nextAxis);
//	//this->greater = new KDBranch(new vector<Triangle>(triangles->begin() + lowerIndex, triangles->end()), nextAxis);
//	/*if (upperIndex < OBJECTS_PER_LEAF) {
//		this->lesser = new KDLeaf(new vector<Triangle>(triangles->begin() + 0, triangles->begin() + upperIndex), nextAxis);
//	}
//	else {
//		this->lesser = new KDBranch(new vector<Triangle>(triangles->begin() + 0, triangles->begin() + lowerIndex), nextAxis, layer + 1);
//	}
//	if (triangles->size() - (lowerIndex + 1) < OBJECTS_PER_LEAF) {
//		this->greater = new KDLeaf(new vector<Triangle>(triangles->begin() + lowerIndex + 1, triangles->end()), nextAxis);
//	}
//	else {
//		this->greater = new KDBranch(new vector<Triangle>(triangles->begin() + lowerIndex, triangles->end()), nextAxis, layer + 1);
//	}*/
//}
//
////TODO: perhaps consider instead of having an intersects, just give those that
////intersect to both left and right, to be sorted at a future time, with only leaf nodes
////having any triangles at all.
//
//KDBranch::KDBranch(vector<Triangle>* triangles, Axis axis = X, int layer = 0) : KDNode()
//{
//	printf("\n======================\nlayer: %i, triangles: %li\n", layer, triangles->size());
//	this->isLeaf = false;
//	//printf("about to hit stoppoint at layer: %i\n", layer);
//
//	this->axis = axis;
//	auto compTri = [axis](const Triangle& a, const Triangle& b)
//	{
//		return int32_t(a.middleBounding[axis] < b.middleBounding[axis]?1:0);
//	};
//
//
//	//printf("here1\n");
//	std::sort(execution::par_unseq, triangles->begin(), triangles->end(), compTri);
//	//printf("here2\n");
//	int32_t midIndex = triangles->size() / 2;
//	//printf("here3\n");
//	
//	this->value = triangles->at(midIndex).middleBounding[axis];
//	int32_t upperIndex = midIndex;
//	int32_t lowerIndex = midIndex;
//	while (upperIndex < triangles->size()-1 && triangles->at(upperIndex).minBounding[axis] < this->value) {
//		upperIndex++;
//	}
//	while (lowerIndex > 0 && triangles->at(lowerIndex).maxBounding[axis] > this->value) {
//		lowerIndex--;
//	}
//
//	Axis nextAxis = Axis((axis + 1) % 3);
//	//this->intersects = vector<Triangle>(triangles.begin() + lowerIndex, triangles.begin() + upperIndex+1);
//
//	printf("upper: %i, lower: %i\n", upperIndex, lowerIndex);
//
//
//
//
//	//TODO: final condition not being reached, need to figure out how to make leaves better
//	if (upperIndex < OBJECTS_PER_LEAF) {
//		printf("making leaf from %i, %i\n", 0, upperIndex);
//		this->lesser = new KDLeaf(new vector<Triangle>(triangles->begin() + 0, triangles->begin() + upperIndex), nextAxis);
//	} else {
//		printf("making branch from %i, %i\n", 0, triangles->end() - (triangles->begin() + lowerIndex));
//		this->lesser = new KDBranch(new vector<Triangle>(triangles->begin() + 0, triangles->begin() + lowerIndex), nextAxis, layer+1);
//	}
//	if (triangles->size()-(lowerIndex+1) < OBJECTS_PER_LEAF) {
//		printf("making lwer leaf from %i, %i\n", lowerIndex + 1, (triangles->end()) - (triangles->begin()));
//		this->greater = new KDLeaf(new vector<Triangle>(triangles->begin() + lowerIndex + 1, triangles->end()), nextAxis);
//	} else {
//		printf("making lower branch from %i, %i\n", lowerIndex, triangles->end() - (triangles->begin()));
//		this->greater = new KDBranch(new vector<Triangle>(triangles->begin() + lowerIndex, triangles->end()), nextAxis, layer+1);
//	}
//}
//
//KDBranch::~KDBranch()
//{
//}
//
//
//enum KDCase {
//	POSITIVE,
//	NEGATIVE,
//	POSTONEG,
//	NEGTOPOS
//};
//
//
//
//KDCase rayCase(const KDBranch* kdn, const Ray& ray) {
//	double rayvalue = ray.origin[kdn->axis];
//	double raydir = ray.direction[kdn->axis];
//	constexpr double epsilon = glm::epsilon<double>();
//	if (rayvalue > kdn->value) {
//		if (raydir <= epsilon) {
//			return POSTONEG;
//		}
//		return POSITIVE;
//	}
//	else {
//		if (raydir >= -epsilon) {
//			return NEGTOPOS;
//		}
//		return NEGATIVE;
//	}
//}
//
//
//
//
//
//
//
////Triangle* evaluateIntersects(const KDNode* kdn, const Ray& ray, vec3& barycentric, double& interDepth) {
//	
//
//	/*vec3 hitPoint;
//	Triangle* closest = NULL;
//	double closestDepth = INFINITY;
//	
//	for (const Triangle& tri : kdn->intersects) {
//		if (rayHit(ray, tri.minBounding, tri.maxBounding)) {
//			if (rayHit(ray, tri, barycentric, hitPoint)) {
//				double depth = glm::length2(hitPoint-ray.origin);
//				if (depth < closestDepth) {
//					closest = (Triangle*) & tri;
//					closestDepth = depth;
//				}
//			}
//		}
//	}
//	interDepth = closestDepth;
//	return closest;*/
//	//return NULL;
////}
//
//
//struct StackEntry {
//	KDNode* node;
//	double tMin;
//	double tMax;
//};
//
//struct KdRayResult {
//	const Triangle* tri;
//	vec3 bary;
//	double depth;
//};
//
///*Based on the paper:
//Interactive k-D Tree GPU Raytracing
//Daniel Reiter Horn Jeremy Sugerman Mike Houston Pat Hanrahan
//Stanford University */
///*KdRayResult kdRayHit(KDNode* kdn, const Ray& ray, vec3& barycentric, double& depthOut, const vec3& sceneMin, const vec3& sceneMax) {
//	stack <StackEntry> kdStack;
//	kdStack.push(StackEntry(kdn, sceneMin.x, sceneMax.x));
//	//double tHit = INFINITY;
//	KdRayResult closest = KdRayResult(NULL, vec3(0), INFINITY);
//	while (!kdStack.empty()) {
//		StackEntry currEntry = kdStack.top();
//		kdStack.pop();
//		while (!currEntry.node->isLeaf) {
//			KDBranch* currBranch = (KDBranch*)currEntry.node;
//			Axis a = currBranch->axis;
//			double tSplit = (currBranch->value - ray.origin[a]) / ray.direction[a];
//			KDNode* first;
//			KDNode* second;
//			if (ray.direction[a] < 0) {
//				first = currBranch->greater;
//				second = currBranch->lesser;
//			} else {
//				first = currBranch->lesser;
//				second = currBranch->greater;
//			}
//			if (tSplit >= currEntry.tMax || tSplit < 0) {
//				currEntry.node = first;
//			} else if (tSplit < currEntry.tMin) {
//				currEntry.node = second;
//			} else {
//				kdStack.push(StackEntry(second, tSplit, currEntry.tMax));
//				currEntry.node = first;
//				currEntry.tMax = tSplit;
//			}
//		}
//		KDLeaf* currLeaf = (KDLeaf*)currEntry.node;
//
//		for (const Triangle& t : currLeaf->tris) {
//			vec3 currBary;
//			double d = kdRayHit(ray, t, currBary);
//			if (d < closest.depth) {
//				closest = KdRayResult(&t, currBary, d);
//			}
//		}
//	}
//	return closest;
//}
//
//
//*/
//
//
//
//
//
//
//
//KdRayResult rayHit(KDNode* kdn, const Ray& ray) {
//	stack <KDNode*> kdStack;
//	kdStack.push(kdn);
//	//double tHit = INFINITY;
//	KdRayResult closest = KdRayResult(NULL, vec3(0), INFINITY);
//	vec3 currBary;
//	while (!kdStack.empty()) {
//		KDNode* currNode = kdStack.top();
//		kdStack.pop();
//		while (currNode != NULL && !currNode->isLeaf) {
//			KDBranch* currBranch = (KDBranch*)currNode;
//			Axis a = currBranch->axis;
//			double tSplit = (currBranch->value - ray.origin[a]) / ray.direction[a];
//			KDNode* first;
//			KDNode* second;
//
//			switch (rayCase(currBranch, ray)) {
//				case POSITIVE:
//					currNode = currBranch->greater;
//					break;
//				case NEGATIVE:
//					currNode = currBranch->lesser;
//					break;
//				case POSTONEG:
//					kdStack.push(currBranch->lesser);
//					currNode = currBranch->greater;
//					break;
//				case NEGTOPOS:
//					kdStack.push(currBranch->greater);
//					currNode = currBranch->lesser;
//					break;
//				default:
//					printf("bad raycase()\n");
//					exit(-1);
//			}
//		}
//		KDLeaf* currLeaf = (KDLeaf*)currNode;
//		if (currLeaf != NULL) {
//			for (const Triangle& t : *(currLeaf->tris)) {
//				double depth = kdRayHit(ray, t, currBary);
//				if (depth < closest.depth) {
//					closest = KdRayResult(&t, currBary, depth);
//				}
//			}
//			if (closest.tri != NULL) {
//				return closest;
//			}
//		}
//	}
//	return closest;
//
//	/*double trueDepth = INFINITY;
//	Triangle* trueT = NULL;
//
//	Triangle* posT = NULL;
//	double posDepth = INFINITY;
//
//	Triangle* interT = NULL;
//	double interDepth = INFINITY;
//
//	Triangle* negT = NULL;
//	double negDepth = INFINITY;*/
//
//	/*if (kdn->greater != NULL) {
//		posT = rayHit(kdn->greater, ray, barycentric, posDepth);
//	}
//	interT = evaluateIntersects(kdn, ray, barycentric, interDepth);
//	if (kdn->lesser != NULL) {
//		negT = rayHit(kdn->lesser, ray, barycentric, negDepth);
//	}*/
//
//	/*if (posDepth < trueDepth && posT != NULL) {
//		trueDepth = posDepth;
//		trueT = posT;
//	}
//	if (negDepth < trueDepth && negT != NULL) {
//		trueDepth = negDepth;
//		trueT = negT;
//	}
//	depthOut = trueDepth;
//	return trueT;*/
//
//
//	/*
//	switch (rayHit(kdn, ray)) {
//		case POSITIVE:
//			if (kdn->greater != NULL) {
//				t = rayHit(kdn->greater, ray, barycentric);
//			}
//			if (t != NULL) {
//				t = evaluateIntersects(kdn, ray, barycentric);
//			}
//			//rayHit() but for the positive side of kdn,
//			//then evaluate the intersect ones for kdn
//			break;
//		case NEGATIVE:
//			if (kdn->lesser != NULL) {
//				t = rayHit(kdn->lesser, ray, barycentric);
//			}
//			if (t == NULL) {
//				t = evaluateIntersects(kdn, ray, barycentric);
//			}
//			//rayHit() but for the negative side of kdn,
//			//then evaluate the intersect ones for kdn
//			break;
//		case POSTONEG:
//			if (kdn->greater != NULL) {
//				t = rayHit(kdn->greater, ray, barycentric);
//			}
//			if (t == NULL) {
//				t = evaluateIntersects(kdn, ray, barycentric);
//			}
//			if (t == NULL && kdn->lesser != NULL) {
//				t = rayHit(kdn->lesser, ray, barycentric);
//			}
//			//rayHit() but for the positive side of kdn,
//			//then evaluate the intersect ones for kdn
//			//rayHit() but for the negative side of kdn,
//			break;
//		case NEGTOPOS:
//			if (kdn->lesser != NULL) {
//				t = rayHit(kdn->lesser, ray, barycentric);
//			}
//			if (t == NULL) {
//				t = evaluateIntersects(kdn, ray, barycentric);
//			}
//			if (t == NULL && kdn->greater != NULL) {
//				t = rayHit(kdn->greater, ray, barycentric);
//			}
//			//rayHit() but for the negative side of kdn,
//			//then evaluate the intersect ones for kdn
//			//rayHit() but for the positive side of kdn,
//			break;
//		default:
//			printf("\n\nBad KDCase\n\n");
//	}
//	return t;*/
//}
//
//
//
//
//
//
//
//#endif;
