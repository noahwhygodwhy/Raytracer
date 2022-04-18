#ifndef KDTREE_H
#define KDTREE_H

#ifndef CPP
#define CPP
#endif

#include <vector>
#include <algorithm>

#include <queue>

#include <glm/glm.hpp>

#include "SharedStructs.cl"

using namespace std;
using namespace glm;


enum Axis {
	X = 0, 
	Y,
	Z
};

//class KDNode
//{
//public:
//	KDNode(bool isLeaf);
//	~KDNode();
//	bool isLeaf;
//private:

//};

class KDNodeNew;


class KDNodeNew
{
public:
	KDNodeNew* greater;
	KDNodeNew* lesser;
	AABB area;
	Axis axis;
	double value;
	vector<int> shapes;

	KDNodeNew(KDNodeNew* greater, KDNodeNew* lesser, AABB area, Axis axis, double value, const vector<int>& shapes);
	KDNodeNew(const vector<int>& shapes);

	~KDNodeNew();

private:

};
//
//class KDLeaf : public KDNode
//{
//public:
//	vector<int> shapes;
//	KDLeaf(vector<int> shapes);
//	~KDLeaf();
//
//private:
//
//};


KDNodeNew* buildKDTree(vector<int> shapeIndexes, vector<UShape> shapes, vector<Vertex> vertices, AABB box, Axis axis = X, int recusionWithoutChange = 0, int layer = 0);


//bool traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTIme, int layer = 0);
//bool rayHitListOfShapes(vector<Shape*> shapes, const Ray& ray, HitResult& hit, double currentTIme);


//vector<GKDNode> transcribeKDTree(KDNode* tree, vector<UShape>& shapes);

#endif