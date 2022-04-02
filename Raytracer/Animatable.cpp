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



dmat4 circle0(double currentTime) {

	currentTime += 0.000;
	return glm::translate(dmat4(1.0), dvec3(sin(currentTime)*15, 15, cos(currentTime)*15));

}dmat4 circle1(double currentTime) {

	currentTime += 0.333*2.0*glm::pi<double>();
	return glm::translate(dmat4(1.0), dvec3(sin(currentTime) * 15, 15, cos(currentTime) * 15));

}dmat4 circle2(double currentTime) {
	currentTime += 0.666 * 2.0 * glm::pi<double>();
	return glm::translate(dmat4(1.0), dvec3(sin(currentTime) * 15, 15, cos(currentTime) * 15));
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
