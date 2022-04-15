//#ifndef KDTREE_H
//#define KDTREE_H
//
//#include <vector>
//#include <algorithm>
//#include "Shape.hpp"
//
//
//
//using namespace std;
//using namespace glm;
//
//
//enum Axis {
//	X = 0, 
//	Y,
//	Z
//};
//
//class KDNode
//{
//public:
//	KDNode(bool isLeaf);
//	~KDNode();
//	bool isLeaf;
//private:
//
//};
//
//
//class KDBranch : public KDNode
//{
//public:
//	KDNode* greater;
//	KDNode* lesser;
//	AABB area;
//	Axis axis;
//	double value;
//	KDBranch(KDNode* greater, KDNode* lesser, AABB area, Axis axis, double value);
//	~KDBranch();
//
//private:
//
//};
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
//
//KDNode* buildKDTree(vector<UShape> shapes, vector<Vertex> vertices, AABB box, Axis axis = X, int recusionWithoutChange = 0, int layer = 0);
////bool traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTIme, int layer = 0);
////bool rayHitListOfShapes(vector<Shape*> shapes, const Ray& ray, HitResult& hit, double currentTIme);
//
//#endif