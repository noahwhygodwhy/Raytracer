//#include "KDTree.hpp"
//
//using namespace std;
//extern bool prd;
//
//
//KDNodeNew::KDNodeNew(KDNodeNew* greater, KDNodeNew* lesser, AABB area, Axis axis, double value, const vector<int>& shapes)
//{
//	this->greater = greater;
//	this->lesser = lesser;
//	this->shapes = shapes;
//	this->area = area;
//	this->axis = axis;
//	this->value = value;
//
//}
//KDNodeNew::KDNodeNew(const vector<int>& shapes) {
//
//	this->greater = greater;
//	this->lesser = lesser;
//	this->shapes = shapes;
//}
//
//
//
//void sortShapesMax(vector<int> shapeIndexes, vector<UShape>& shapes, Axis axis) {
//	sort(shapeIndexes.begin(), 
//		shapeIndexes.end(),
//		[axis, shapes](int a, int b) -> bool {
//
//			return shapes.at(b).boundingBox.max[axis] > shapes.at(a).boundingBox.max[axis]; 
//		});
//}
//void sortShapesMin(vector<int> shapeIndexes, vector<UShape>& shapes, Axis axis) {
//	sort(shapeIndexes.begin(),
//		shapeIndexes.end(),
//		[axis, shapes](int a, int b) -> bool {
//
//			return shapes.at(b).boundingBox.min[axis] > shapes.at(a).boundingBox.min[axis];
//		});
//}
//
//
//
//
//
//const bool dKD = false;
//
//KDNodeNew* buildKDTree(vector<int> shapeIndexes, vector<UShape> shapes, vector<Vertex> vertices, AABB box, Axis axis, int recursionWithoutChange, int layer) {
//	if (shapeIndexes.size() < 7 || recursionWithoutChange > 3 || layer > 20) {//TODO: this is basic, and can be optimized
//		return new KDNodeNew(shapeIndexes);
//	}
//	AABB midBounding = shapes.at(shapeIndexes.at(shapeIndexes.size() / 2u)).boundingBox;
//	double value = (midBounding.min[axis] + midBounding.max[axis]) / 2.0f;
//	Axis nextAxis = Axis((axis + 1) % 3);
//
//
//
//	int maxLesserIndex = 0;
//	sortShapesMin(shapeIndexes, shapes, axis);
//
//	while (maxLesserIndex < shapeIndexes.size() && shapes.at(shapeIndexes.at(maxLesserIndex)).boundingBox.max[axis] <= value) {
//		maxLesserIndex++;
//	}
//	KDNodeNew* lesserNode = NULL;
//	if (maxLesserIndex > 0) {
//		vector<int> lesserShapes(shapeIndexes.begin(), shapeIndexes.begin() + maxLesserIndex);
//		AABB lesserBox = box;
//		lesserBox.max[axis] = value;
//		int newRWC = 0;
//		if (lesserShapes.size() == shapes.size()) {
//			newRWC = recursionWithoutChange + 1;
//		}
//		lesserNode = buildKDTree(lesserShapes, shapes, vertices, lesserBox, nextAxis, newRWC, layer + 1);
//	}
//	else {
//		lesserNode = NULL;
//	}
//
//
//
//	int minGreaterIndex = shapeIndexes.size() - 1;
//	sortShapesMax(shapeIndexes, shapes, axis);
//	while (minGreaterIndex > 0 && shapes.at(shapeIndexes.at(minGreaterIndex)).boundingBox.min[axis] >= value ) {
//		minGreaterIndex--;
//	}
//	KDNodeNew* greaterNode = NULL;
//	if (minGreaterIndex < (shapeIndexes.size() - 1)) {
//		vector<int> greaterShapes(shapeIndexes.begin() + minGreaterIndex, shapeIndexes.end());
//		AABB greaterBox = box;
//		greaterBox.min[axis] = value;
//		int newRWC = 0;
//		if (greaterShapes.size() == shapes.size()) {
//			newRWC = recursionWithoutChange + 1;
//		}
//		greaterNode = buildKDTree(greaterShapes, shapes, vertices, greaterBox, nextAxis, newRWC, layer + 1);
//	}
//	else {
//		greaterNode = NULL;
//	}
//	//TODO:
//	vector<int> middleShapes(shapeIndexes.begin() + minGreaterIndex + 1, shapeIndexes.begin() + maxLesserIndex - 1);
//
//	return new KDNodeNew(
//		greaterNode,
//		lesserNode,
//		box,
//		axis,
//		value,
//		middleShapes
//	);
//}
//
//vector<GKDNode> transcribeKDTree(KDNodeNew* tree, vector<UShape>& shapes) {
//	queue<KDNodeNew*> nodeQ;
//	nodeQ.push(tree);
//	int64_t i = 0;
//	while (!nodeQ.empty()) {
//
//		KDNodeNew* currNode = nodeQ.front();
//		nodeQ.pop();
//		if(currNode->lesser != 0){
//			nodeQ.push(currNode->lesser);
//		}
//		if (currNode->greater != 0) {
//			nodeQ.push(currNode->greater);
//		}
//
//
//		/*typedef struct alignas(16) GKDNode
//		{
//			AABB area;
//			float value;
//			int greater;
//			int lesser;
//			int axis;
//		} GKDNode;*/
//
//
//
//
//
//	}
//
//}
//
//
//
////https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.252.9650&rep=rep1&type=pdf
////bool traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTime, int layer) {
////	if (tree->isLeaf) {
////		bool itHit =  rayHitListOfShapes(
////			((KDLeaf*)tree)->shapes,
////			ray, hit, currentTime);
////		return itHit;
////	}
////	else {
////		KDBranch* currBranch = (KDBranch*)tree;
////		optional<dvec3> enter;
////		optional<dvec3> exit;
////		bool didItHit = currBranch->area.rayAABB(ray, enter, exit);
////		if (didItHit) {
////			if (enter.value()[currBranch->axis] <= currBranch->value) {
////
////				if (exit.value()[currBranch->axis] < currBranch->value) {
////					return traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
////				}
////				else {
////					bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
////					bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
////					return leftSideHit || rightSideHit;
////					/*bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
////					/*if (leftSideHit) {
////						return leftSideHit;
////					}
////					return traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);*/
////				}
////			}
////			else {
////				if (exit.value()[currBranch->axis] > currBranch->value) {
////					return traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
////				}
////				else {
////					bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
////					bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
////					return leftSideHit || rightSideHit;
////					/*bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
////					if (rightSideHit) {
////						return rightSideHit;
////					}
////					return traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);*/
////				}
////			}
////		}
////	}
////}
//// 