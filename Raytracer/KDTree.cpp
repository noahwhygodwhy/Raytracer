#include "KDTree.hpp"
#include "SharedStructs.cl"
extern bool prd;

KDNode::KDNode(bool isLeaf)
{
	this->isLeaf = isLeaf;

}
KDBranch::KDBranch(KDNode* greater, KDNode* lesser, AABB area, Axis axis, double value) : KDNode(false)
{
	this->greater = greater;
	this->lesser = lesser;
	this->area = area;
	this->axis = axis;
	this->value = value;

}
KDLeaf::KDLeaf(vector<int> shapes):KDNode(true)
{
	this->shapes = shapes;
}

KDNode::~KDNode()
{
}

KDBranch::~KDBranch()
{
}

KDLeaf::~KDLeaf()
{
}



/*double compX(Shape* a, Shape* b) {
	return a->boundingBox.mid - b->boundingBox.mid;
}
double compShape(Shape* a, Shape* b) {
	return a->boundingBox.mid - b->boundingBox.mid;
}*/

void sortShapesMax (vector<Shape*>& shapes, Axis axis) {
	sort(shapes.begin(), 
		shapes.end(), 
		[axis](Shape* a, Shape* b) -> bool {

			return b->boundingBox.max[axis] > a->boundingBox.max[axis]; 
		});
}
void sortShapesMin(vector<Shape*>& shapes, Axis axis) {
	sort(shapes.begin(),
		shapes.end(),
		[axis](Shape* a, Shape* b) -> bool {

			return b->boundingBox.min[axis] > a->boundingBox.min[axis];
		});
}





const bool dKD = false;

KDNode* buildKDTree(vector<Shape*> shapes, vector<Vertex> vertices, AABB box, Axis axis, int recursionWithoutChange, int layer) {
	if (shapes.size() < 7 || recursionWithoutChange > 3 || layer > 20) {//TODO: this is basic, and can be optimized
		return new KDLeaf(shapes); 
	}
	double value = shapes.at(shapes.size() / 2u)->boundingBox.mid[axis];
	Axis nextAxis = Axis((axis + 1) % 3);



	int maxLesserIndex = 0;
	sortShapesMin(shapes, axis);
	while (maxLesserIndex < shapes.size() && shapes.at(maxLesserIndex)->boundingBox.min[axis] <= value) {
		maxLesserIndex++;
	}
	KDNode* lesserNode = NULL;
	if (maxLesserIndex > 0) {
		vector<Shape*> lesserShapes(shapes.begin(), shapes.begin() + maxLesserIndex);
		AABB lesserBox = box;
		lesserBox.max[axis] = value;
		int newRWC = 0;
		if (lesserShapes.size() == shapes.size()) {
			newRWC = recursionWithoutChange + 1;
		}
		lesserNode = buildKDTree(lesserShapes, lesserBox, nextAxis, newRWC, layer + 1);
	}
	else {
		lesserNode = new KDLeaf({});
	}



	int minGreaterIndex = shapes.size() - 1;
	sortShapesMax(shapes, axis);
	while (minGreaterIndex > 0 && shapes.at(minGreaterIndex)->boundingBox.max[axis] >= value ) {
		minGreaterIndex--;
	}
	KDNode* greaterNode = NULL;
	if (minGreaterIndex < (shapes.size() - 1)) {
		vector<Shape*> greaterShapes(shapes.begin() + minGreaterIndex, shapes.end());
		AABB greaterBox = box;
		greaterBox.min[axis] = value;
		int newRWC = 0;
		if (greaterShapes.size() == shapes.size()) {
			newRWC = recursionWithoutChange + 1;
		}
		greaterNode = buildKDTree(greaterShapes, greaterBox, nextAxis, newRWC, layer + 1);
	}
	else {
		greaterNode = new KDLeaf({});
	}

	return new KDBranch(
		greaterNode,
		lesserNode,
		box,
		axis,
		value
	);
}



//https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.252.9650&rep=rep1&type=pdf
bool traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTime, int layer) {
	if (tree->isLeaf) {
		bool itHit =  rayHitListOfShapes(
			((KDLeaf*)tree)->shapes,
			ray, hit, currentTime);
		return itHit;
	}
	else {
		KDBranch* currBranch = (KDBranch*)tree;
		optional<dvec3> enter;
		optional<dvec3> exit;
		bool didItHit = currBranch->area.rayAABB(ray, enter, exit);
		if (didItHit) {
			if (enter.value()[currBranch->axis] <= currBranch->value) {

				if (exit.value()[currBranch->axis] < currBranch->value) {
					return traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
				}
				else {
					bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
					bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
					return leftSideHit || rightSideHit;
					/*bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
					/*if (leftSideHit) {
						return leftSideHit;
					}
					return traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);*/
				}
			}
			else {
				if (exit.value()[currBranch->axis] > currBranch->value) {
					return traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
				}
				else {
					bool leftSideHit = traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);
					bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
					return leftSideHit || rightSideHit;
					/*bool rightSideHit = traverseKDTree(currBranch->greater, ray, hit, currentTime, layer + 1);
					if (rightSideHit) {
						return rightSideHit;
					}
					return traverseKDTree(currBranch->lesser, ray, hit, currentTime, layer + 1);*/
				}
			}
		}
	}
}
 