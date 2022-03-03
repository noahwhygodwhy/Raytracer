#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <algorithm>
#include "Shape.hpp"


using namespace std;
using namespace glm;


enum Axis {
	X = 0, 
	Y,
	Z
};

class KDNode
{
public:
	KDNode();
	~KDNode();
	bool isLeaf;
private:

};


class KDBranch : public KDNode
{
public:
	KDNode* greater;
	KDNode* lesser;
	AABB area;
	Axis axis;
	double value;
	KDBranch(KDNode* greater, KDNode* lesser, AABB area, Axis axis, double value);
	~KDBranch();

private:

};

class KDLeaf : public KDNode
{
public:
	vector<Shape*> shapes;
	KDLeaf(vector<Shape*> shapes);
	~KDLeaf();

private:

};

KDNode* buildKDTree(vector<Shape*> shapes, AABB box, Axis axis = X);
void traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTIme);
void rayHitListOfShapes(vector<Shape*> shapes, const Ray& ray, HitResult& hit, double currentTIme);

#endif