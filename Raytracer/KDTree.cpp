#include "KDTree.hpp"

KDNode::KDNode()
{

}
KDBranch::KDBranch(KDNode* greater, KDNode* lesser, AABB area, Axis axis, double value)
{
	this->greater = greater;
	this->lesser = lesser;
	this->area = area;
	this->axis = axis;
	this->value = value;
	this->isLeaf = false;

}
KDLeaf::KDLeaf(vector<Shape*> shapes)
{
	this->shapes = shapes;
	this->isLeaf = true;
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

void sortShapes (vector<Shape*>& shapes, Axis axis) {
	sort(shapes.begin(), 
		shapes.end(), 
		[axis](Shape* a, Shape* b) -> bool {
			return a->boundingBox.mid[axis] < b->boundingBox.mid[axis]; 
		});
}


//https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.252.9650&rep=rep1&type=pdf
KDNode* buildKDTree(vector<Shape*> shapes, AABB box, Axis axis) {
	if (shapes.size() < 50) {
		return new KDLeaf(shapes);
	}
	sortShapes(shapes, axis);
	double value = shapes.at(shapes.size() / 2)->boundingBox.mid[axis];
	size_t maxLesserIndex = 0;
	size_t minGreaterIndex = shapes.size() - 1;
	while (shapes.at(maxLesserIndex)->boundingBox.max[axis] >= value) {
		maxLesserIndex++;
	}
	while (shapes.at(minGreaterIndex)->boundingBox.min[axis] >= value) {
		minGreaterIndex--;
	}
	vector<Shape*> greaterShapes(shapes.begin()+minGreaterIndex, shapes.end());
	vector<Shape*> lesserShapes(shapes.begin(), shapes.begin() + maxLesserIndex);
	AABB greaterBox = box;
	greaterBox.min[axis] = value;
	AABB lesserBox = box;
	greaterBox.max[axis] = value;

	Axis nextAxis = Axis((axis + 1) % 3);

	return new KDBranch(
		buildKDTree(greaterShapes, greaterBox, nextAxis),
		buildKDTree(lesserShapes, lesserBox, nextAxis),
		box,
		axis,
		value
	);

}

void traverseKDTree(KDNode* tree, const Ray& ray, HitResult& hit, double currentTime) {
	if (tree->isLeaf) {
		rayHitListOfShapes(
			((KDLeaf*)tree)->shapes,
			ray, hit, currentTime);
	}
	else {
		
	}
}

void rayHitListOfShapes(vector<Shape*> shapes, const Ray& ray, HitResult& minRayResult, double currentTime) {
	for (Shape* shape : shapes) {
		HitResult rayResult;
		if (shape->rayAABB(ray)) {
			if (shape->rayHit(ray, rayResult, currentTime)) {
				if (rayResult.depth < minRayResult.depth) {
					minRayResult = rayResult;
				}
			}
		}
	}
}