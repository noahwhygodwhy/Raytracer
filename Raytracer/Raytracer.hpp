#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <iostream>
#include <vector>
#include <fstream>
#include <mutex>
#include <execution>
#include <ppl.h>
#include <thread>
#include <filesystem>
#include <ctime>
#include <random>
//#include <queue>
//#include <execution>

#include <CL/cl.h>


#include "glad/glad.h"
#include <GLFW/glfw3.h>
#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image_write.h"




#include "Shader.hpp"
#include "Model.hpp"
//#include "WorkerPool.hpp"

#include "CoordinateHelpers.hpp"

#include "Material.hpp"
#include "Texture.hpp"
#include "Procedural.hpp"

#include "Triangle.hpp"
#include "Ray.hpp"
#include "KDTree.hpp"
#include "Shape.hpp"
#include "Sphere.hpp"
#include "Biconvex.hpp"



using namespace std;
using namespace std::filesystem;
using namespace glm;


struct FrameInfo {
	vector<Shape*> shapes;
	//vector<Light*> lights;
	KDNode* kdTree;
	//dmat4 view;
	dvec3 camPosition;
	double currentTime;
};

void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

}


fvec3* frameBuffer;
fvec3* drawBuffer;


#endif