#include "KDTree.hpp"

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
KDLeaf::KDLeaf(vector<Shape*> shapes):KDNode(true)
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

KDNode* buildKDTree(vector<Shape*> shapes, AABB box, Axis axis, int recursionWithoutChange, int layer) {
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


void printKDBranch(KDBranch* b) {
	printf("%s split at %f\n", b->axis == 0 ? "X" : b->axis == 1 ? "Y" : "Z", b->value);
}
void spaces(int layer) {
	for (int i = 0; i < layer; i++) {
		printf(" ");
	}
}
void printKDTree(KDNode* b, int layer) {
	//spaces(layer);
	if (b->isLeaf) {
		KDLeaf* l = (KDLeaf*)b;
		printf("leaf with %zu shapes\n", l->shapes.size());
	}
	else {
		KDBranch* l = (KDBranch*)b;

		printf("branch axis %s split at %f\n", l->axis == 0 ? "X" : l->axis == 1 ? "Y" : "Z", l->value);
		spaces(layer+1);
		printf("lesser: ");
		printKDTree(l->lesser, layer + 1);
		spaces(layer+1);
		printf("greater: ");
		printKDTree(l->greater, layer + 1);
		
	}
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

bool rayHitListOfShapes(vector<Shape*> shapes, const Ray& ray, HitResult& minRayResult, double currentTime) {
	bool toReturn = false;
	if (prd)printf("rayhitlistofshapes with ray: %s, %s\n", glm::to_string(ray.origin).c_str(), glm::to_string(ray.direction).c_str());
	if(prd)printf("shapes: %zu\n", shapes.size());
	for (Shape* shape : shapes) {
		HitResult rayResult;
		if (shape->rayAABB(ray)) {//TODO reenable this
			if (prd)printf("rayaabb\n");
			if (shape->rayHit(ray, rayResult, currentTime)) {
				if(prd)printf("ray hit at %s\n", glm::to_string(rayResult.position).c_str());
				if(prd)printf("depth of %f\n", rayResult.depth);
				if (rayResult.depth < minRayResult.depth) {
					toReturn = true;
					minRayResult = rayResult;
				}
			}
		}
	}
	return toReturn;
}