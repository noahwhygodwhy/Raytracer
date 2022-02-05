#ifndef ANIMATABLE_H
#define ANIMATABLE_H

#include <glm/glm.hpp>
using namespace std;
using namespace glm;


dmat4 noMovement(double currentTime) {
	return dmat4(1.0);
}

dmat4 oscilateX(double currentTime) {
	return glm::translate(dmat4(1.0), dvec3(sin(currentTime) * 5.0, 0.0, 0.0));
}



class Animatable
{
public:
	Animatable(dmat4(*movementFunction)(double currentTime));
	~Animatable();
	dmat4 getModel(double currentTime) const;
private:
	dmat4 (*movementFunction)(double currentTime);//double being the current frame time in seconds

};

Animatable::Animatable(dmat4(*movementFunction)(double currentTime) = noMovement)
{
	this->movementFunction = movementFunction;
}

Animatable::~Animatable()
{
}

dmat4 Animatable::getModel(double currentTime) const {
	return this->movementFunction(currentTime);
}





#endif