#ifndef ANIMATABLE_H
#define ANIMATABLE_H


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace std;
using namespace glm;


dmat4 noMovement(double currentTime);

dmat4 oscilateX(double currentTime);
dmat4 oscilateY(double currentTime);
dmat4 rotateY(double currentTime);

class Animatable
{
public:
	Animatable(dmat4(*movementFunction)(double currentTime));
	~Animatable();
	dmat4 getModel(double currentTime) const;
private:
	dmat4 (*movementFunction)(double currentTime);//double being the current frame time in seconds

};





#endif