#include "Animatable.hpp"
using namespace std;
using namespace glm;

dmat4 noMovement(double currentTime) {
	return dmat4(1.0);
}

dmat4 oscilateX(double currentTime) {
	return glm::translate(dmat4(1.0), dvec3(sin(currentTime) * 5.0, 0.0, 0.0));
}
dmat4 oscilateY(double currentTime) {
	return glm::translate(dmat4(1.0), dvec3(0.0, sin(currentTime) * 5.0, 0.0));
}
dmat4 rotateY(double currentTime) {
	return glm::rotate(dmat4(1.0), currentTime / 3.0, dvec3(0.0, 1.0, 0.0));

}
Animatable::Animatable(dmat4(*mf)(double currentTime))
{
	this->movementFunction = mf;
}

Animatable::~Animatable()
{
}

dmat4 Animatable::getModel(double currentTime) const {
	return this->movementFunction(currentTime);
}
