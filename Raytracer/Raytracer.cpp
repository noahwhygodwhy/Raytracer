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

//#include <CL/cl.h>


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

#include "Light.hpp"
#include "PointLight.hpp"
#include "DirectionalLight.hpp"
#include "SquareLight.hpp"


using namespace std;
using namespace std::filesystem;
using namespace glm;


#define MAX_PATH 50
#define MAXBOUNCES 3u
//#define OUTPUTFRAMES 151
//#define EVERYFRAME INFINITY
#define CONCURRENT_FOR
#define KDTRACE
#define CIN
#define LIGHTING
#define GLOBAL_LIGHTING

#define LIGHT_MULTISAMPLES_N 1
#define PIXEL_MULTISAMPLE_N 1
#define MONTE_CARLO_SAMPLES 15000

#define PATH_TRACE


#define LIGHT_CONTSTANT 1.0
#define LIGHT_LINEAR 0.09
#define LIGHT_QUADRATIC 0.032
double getAttenuation(double distance) {
	return (LIGHT_CONTSTANT + (LIGHT_LINEAR * distance) + (LIGHT_QUADRATIC * (distance * distance)));
}


bool prd = false; //print debuging for refraction

uint32_t frameX = 500;
uint32_t frameY = 500;
double frameRatio = double(frameX) / double(frameY);


//dvec3 clearColor(0.21, 0.78, 0.95);
dvec3 clearColor(0.0, 0.0, 0.0);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame
string saveFileDirectory = "";


constexpr double bias = 1e-3;


struct FrameInfo {
	vector<Shape*> shapes;
	vector<Light*> lights;
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

vec3* frameBuffer;
vec3* drawBuffer;

void frameBufferSizeCallback(GLFWwindow* window, uint64_t width, uint64_t height) {
	glViewport(0, 0, GLsizei(width), GLsizei(height));
}

void clearBuffers() {
	for (uint64_t i = 0; i < frameX * frameY; i++) {
		frameBuffer[i] = clearColor;
	}
}
void clearDrawBuffer() {
	for (uint64_t i = 0; i < frameX * frameY; i++) {
		drawBuffer[i] = clearColor;
	}
}

bool frontFacing(vec3 a, vec3 b, vec3 c) {
	mat3 m(a.x, b.x, c.x, a.y, b.y, c.y, 1.0f, 1.0f, 1.0f);
	return glm::determinant(m) > 0;
}

//random double from 0 to 1
double randDubTwo() {
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
dvec3 randomHemisphericalVector(dvec3 normal) {

	normal = glm::normalize(normal);
	dvec3 worldUp = normal.y > 0.9 ? dvec3(0.0, 1.0, 0.0) : dvec3(1.0, 0.0, 0.0);

	double theta = 2 * glm::pi<double>() * randDubTwo();
	double phi = glm::acos(1.0 - (2 * randDubTwo()));

	dvec3 b1 = glm::normalize(glm::cross(normal, worldUp));
	dvec3 b2 = glm::cross(b1, normal);

	dvec3 newVec = dvec3(sin(theta) * sin(phi), sin(theta) * cos(phi), sin(phi));


	dvec3 rotVector = glm::cross(normal, worldUp);
	double rotRads = glm::acos(glm::dot(normal, worldUp));
	return glm::rotate(newVec, rotRads, rotVector);
}

//since it's not part of the assignment to write a function to save to image, i just copied this directly from
//https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/
void saveImage(string filepath, GLFWwindow* w) {

	string outDir = "out/"+saveFileDirectory+"/";
	int width, height;
	glfwGetFramebufferSize(w, &width, &height);
	GLsizei nrChannels = 3;
	GLsizei stride = nrChannels * width;
	stride += (stride % 4) ? (4 - stride % 4) : 0;
	GLsizei bufferSize = stride * height;
	std::vector<char> buffer(bufferSize);
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
	stbi_flip_vertically_on_write(true);
	stbi_write_png((outDir + filepath).c_str(), width, height, nrChannels, buffer.data(), stride);
}




double shlicksApprox(double matIOR, vec3 normal, const vec3& viewVector, bool entering) {

	double prevIOR = 1.0;
	double newIOR = matIOR;
	if (!entering) {
		normal = -normal;
		swap(prevIOR, newIOR);
	}

	double r0 = glm::pow((prevIOR - newIOR) / (prevIOR + newIOR), 2.0);

	double rTheta = r0 + ((1 - r0) * (1 - glm::acos(glm::dot(normal, viewVector))));
	return rTheta;
}

HitResult shootRay(const Ray& ray, const FrameInfo& fi) {

	HitResult minRayResult;
	minRayResult.shape = NULL;
	minRayResult.depth = INFINITY;

#ifdef KDTRACE
	traverseKDTree(fi.kdTree, ray, minRayResult, fi.currentTime);
#else
	rayHitListOfShapes(fi.shapes, ray, minRayResult, fi.currentTime);
#endif
	return minRayResult;
}

dvec3 getRefractionRay(dvec3 hitNormal, dvec3 incidentVector, double objectIOR, bool entering, bool& internalOnly) {


	//TODO: this is where you reverse the normal if the dot of normal and incident vector is above a certain amount
	double closeness = glm::dot(hitNormal, incidentVector);
	double prevIOR = 1.0;
	double newIOR = objectIOR;


	if (!entering) {
		hitNormal = -hitNormal;
		swap(prevIOR, newIOR);
	}

	//double cosA1 = glm::clamp(glm::dot(incidentVector, hitNormal), -1.0, 1.0);
	double cosA1 = glm::dot(incidentVector, hitNormal);

	double sinA1 = glm::sqrt(1.0 - (cosA1 * cosA1));

	double IORRatio = prevIOR / newIOR;

	double sinA2 = sinA1 * IORRatio;

	dvec3 trueReflectDir = incidentVector;
	if (sinA2 <= -1.0 || sinA2 >= 1.0) {//TODO this isn't handled correctly
		internalOnly = true;
		return trueReflectDir;
	}

	double maxCloseness = -INFINITY;
	double k1 = NAN;
	double k2 = NAN;
	solveQuadratic(1.0, 2.0*cosA1, 1.0-(1.0/(IORRatio * IORRatio)), k1, k2);
	// for the solution that k that causes the new ray
	//that has the smallest angle difference between it and the incident ray
	if (!isnan(k1)) {
		dvec3 reflectDir1 = glm::normalize(incidentVector + (k1 * hitNormal));
		double closeness1 = glm::dot(incidentVector, reflectDir1);
		if (closeness1 > maxCloseness && closeness1>=0.0) {
			maxCloseness = closeness1;
			trueReflectDir = reflectDir1;
		}
	}
	if (!isnan(k2)) {
		dvec3 reflectDir2 = glm::normalize(incidentVector + (k2 * hitNormal));
		double closeness2 = glm::dot(incidentVector, reflectDir2);
		if (closeness2 > maxCloseness && closeness2>=0.0) {
			maxCloseness = closeness2;
			trueReflectDir = reflectDir2;
		}
	}

	if (maxCloseness <= 0.0) {
		printf("error calculating refraction angle\n");
		exit(-1);
	}

	double cosA2 = glm::sqrt(1.0 - (sinA2 * sinA2));

	if (cosA1 < 0.0) {
		cosA2 *= -1.0;
	}

	return glm::normalize(trueReflectDir);
}






constexpr double kj = 0.1;//given in the paper
constexpr double kf = 1.12;
constexpr double kg = 1.01;

double F(double x) {
	float num = (1 / (glm::pow(x - kf, 2.0))) - (1 / glm::pow(kf, 2.0));
	float denom = (1 / glm::pow(1 - kf, 2.0)) - (1 / glm::pow(kf, 2.0));
	return num / denom;
}

double G(double x) {
	double num = (1 / glm::pow(1 - kg, 2.0)) - (1 / glm::pow(x - kg, 2.0));
	double denom = (1 / glm::pow(1 - kg, 2.0)) - (1 / glm::pow(kg, 2.0));
	return num / denom;
}




dvec3 lightingFunction(const Ray& originalRay, const Ray& lightRay, const HitResult& minRayResult, const double attenuation, const Material& mat, const dvec3& lightColor) {


	dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal) - lightRay.direction);
	dvec3 H = glm::normalize(lightRay.direction + originalRay.inverseDirection);

	double spec = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, originalRay.inverseDirection)), mat.getNS(minRayResult.uv));//to the specular exponent
	dvec3 specular = (lightColor * spec) / attenuation;
	double diff = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
	dvec3 diffuse = (mat.getColor(minRayResult.uv) * lightColor * diff) / attenuation;
	return (diffuse + specular);
}

dvec3 randomHemisphereVector(double u1, double u2)
{
	double r = glm::sqrt(1.0f - u1 * u1);
	double phi = 2 * glm::pi<double>() * u2;

	return dvec3(cos(phi) * r, sin(phi) * r, u1);
}



dvec3 pathTrace(const Ray& ray, const FrameInfo& fi, double currentIOR = 1.0, uint32_t layer = 0) {

	if (layer > MAX_PATH) {
		if (prd)printf("too deep of a path\n");
		return clearColor * 0.1;
	}
	HitResult minRayResult = shootRay(ray, fi);

	if (minRayResult.shape == NULL) {
		if (prd)printf("didn't hit nuffin\n");
		return clearColor * 0.1;
	}

	dvec2 uv = minRayResult.uv;
	Material mat = minRayResult.shape->mat;

	double trans = mat.getTransparency(uv);
	double smooth = mat.getSmoothness(uv);

	double transparencyDecider = randDubTwo();
	double reflectanceDecider = randDubTwo();


	dvec3 newRayDirection;
	dvec3 downstreamRadiance;


	double hitAngle = glm::acos(glm::dot(minRayResult.normal, ray.inverseDirection));
	bool entering = hitAngle < (glm::pi<double>() / 2.0);
	bool internalOnly;

	dvec3 newRayDir;
	dvec3 newRayPos;


	if (transparencyDecider < trans) {
		if (prd)printf("case1\n");

		newRayDir = getRefractionRay(glm::normalize(minRayResult.normal), glm::normalize(ray.direction), mat.getNI(minRayResult.uv), entering, internalOnly);
		newRayPos = minRayResult.position + (minRayResult.normal * (entering ? -1.0 : 1.0) * bias);
	}
	else if (reflectanceDecider < smooth) {
		if (prd)printf("case2\n");
		newRayDir = -rotate(ray.direction, glm::pi<double>(), minRayResult.normal);
		newRayPos = minRayResult.position + (minRayResult.normal * bias);
	} else {
		if (prd)printf("case3\n");
		//do {
	//		newRayDir = glm::normalize(dvec3(randDubTwo(), randDubTwo(), randDubTwo()));
		//} while (glm::dot(newRayDir, minRayResult.normal) < 0);

		newRayDir = randomHemisphereVector(randDubTwo(), randDubTwo());

		newRayPos = minRayResult.position + (minRayResult.normal * bias);
	}
	//newRayDir = glm::normalize(dvec3(-1.0, 0.005, 0.0));
	//if (prd)printf("reflectionRay: %s, %s\n", glm::to_string(newRayPos).c_str(), glm::to_string(newRayDir).c_str());
	Ray reflectionRay(newRayPos, newRayDir);


	dvec3 thisRadiance = mat.getEmission();
	if (thisRadiance == dvec3(0.0, 0.0, 0.0)) {
		if(prd)printf("radiance is 0\n");
		
		downstreamRadiance = pathTrace(reflectionRay, fi, mat.getNI(uv), layer + 1);

		//calculate this objects emmision, both a reflection of the material, and any emmision

		//TODO: this needs to be an actual lighting BDRF, not just this

//blinn phong (ish)
		dvec3 R = glm::rotate(newRayDir, glm::pi<double>(), minRayResult.normal);
		dvec3 V = ray.origin;
		dvec3 N = minRayResult.normal;
		dvec3 L = newRayDir;

		dvec3 H = (L + V) / glm::length(L + V);

		double diff = glm::dot(L, N);
		double spec = glm::dot(N, H);

		double ks = smooth;
		double kd = 1.0 - ks;

		dvec3 diffuse = kd * diff * downstreamRadiance * mat.getColor(uv);
		dvec3 specular = ks * spec * downstreamRadiance * mat.getColor(uv);

		thisRadiance += (downstreamRadiance* mat.getColor(uv));

	}

	if (prd)printf("\n======================\npath tracing layer %u\n", layer);

	if (prd)printf("transD: %f, reflect: %f\n", transparencyDecider, reflectanceDecider);
	if (prd)printf("reflectionRay: %s, %s\n", glm::to_string(newRayPos).c_str(), glm::to_string(newRayDir).c_str());

	if (prd)printf("thisradiance: %s\n", glm::to_string(thisRadiance).c_str());

	if (prd)printf("returning %s\n", glm::to_string(thisRadiance).c_str());
	if(prd)printf("min ray result depth: %f\n", minRayResult.depth);
	//thisRadiance = thisRadiance * (1.0 / (minRayResult.depth * glm::max(minRayResult.depth, 1.0)));
	if (prd)printf("no instead %s\n", glm::to_string(thisRadiance).c_str());

	return thisRadiance;
	



	/*
		shoot the ray
		then you have the emmision of what it hit, and the distance over which that emmision travels
		calculate the effect of that returned light on this object
		calculate the local lighting situation of this object

		rturn those two as the radiance result
			
	*/

}

bool rayTrace(const Ray& ray, const FrameInfo& fi, dvec3& colorResult, HitResult& minRayResult, double currentIOR = 1.0, uint32_t layer = 0)
{


	if (layer > MAXBOUNCES) {
		//printf("terminating on max bounces\n");
		colorResult = clearColor;
		return false;
	}

	colorResult = dvec3(0);
	minRayResult = shootRay(ray, fi); //ray is in camera space
	if (minRayResult.shape != NULL) {
		
		//if(prd)printf("hit position: %s\n", glm::to_string(minRayResult.position).c_str(), glm::to_string(minRayResult.normal).c_str())
		
		Material mat = minRayResult.shape->mat;

		double hitAngle = glm::acos(glm::dot(minRayResult.normal, ray.inverseDirection));

		bool entering = hitAngle < (glm::pi<double>() / 2.0);


		dvec3 directLightColorResult = dvec3(0.0);
		dvec3 globalLightColorResult = dvec3(0.0);

#ifdef LIGHTING
		for (Light* light : fi.lights) {
#ifdef LIGHT_MULTISAMPLES_N
			int samples = light->multisample ? LIGHT_MULTISAMPLES_N *LIGHT_MULTISAMPLES_N : 1;
#else
			int samples = 1;
#endif

			for (int i = 0; i < samples; i++) {
				Ray lightRay = light->getRay(minRayResult.position + (minRayResult.normal * bias));//should be in camera space
				//if(prd) printf("light ray origin: %s, ray: %s\n", glm::to_string(lightRay.origin).c_str(), glm::to_string(lightRay.direction).c_str());
				HitResult lightRayResult = shootRay(lightRay, fi);
				double lightDistance = light->getDistance(lightRay);

				if (lightDistance < lightRayResult.depth) {//shadows
					double attenuation = light->getAttenuation(lightDistance, LIGHT_CONTSTANT, LIGHT_LINEAR, LIGHT_QUADRATIC);
					directLightColorResult += lightingFunction(ray, lightRay, minRayResult, attenuation, mat, light->color);
				}
			}

			directLightColorResult /= double(samples);
#ifdef GLOBAL_LIGHTING
			int globalSamples = MONTE_CARLO_SAMPLES * ((double(MAXBOUNCES)/double(layer+1))/double(MAXBOUNCES));
			if (prd)printf("layer %u, global samples: %i\n", layer, globalSamples);
			for (int i = 0; i < globalSamples; i++) {
				dvec3 origin = minRayResult.position;
				dvec3 rayDirection = randomHemisphericalVector(minRayResult.normal);
				Ray lightRay(origin, rayDirection);

				if(prd) printf("ray: origin: %s, direction: %s\n", glm::to_string(origin).c_str(), glm::to_string(rayDirection).c_str());
				HitResult lightRayResult;
				dvec3 colorResult;
				bool lightRayHit = rayTrace(lightRay, fi, colorResult, lightRayResult, mat.getNI(minRayResult.uv), layer + 1);

				double rayDistance = glm::distance(origin, lightRayResult.position);
				double attenuation = light->getAttenuation(rayDistance, LIGHT_CONTSTANT, LIGHT_LINEAR, LIGHT_QUADRATIC);
				if (lightRayHit) {
					double attenuation = light->getAttenuation(rayDistance, LIGHT_CONTSTANT, LIGHT_LINEAR, LIGHT_QUADRATIC);

					globalLightColorResult += lightingFunction(ray, lightRay, minRayResult, attenuation, mat, light->color);
				}
				else { 
					globalLightColorResult += clearColor * 0.04;
				}
			}
			globalLightColorResult /= double(MONTE_CARLO_SAMPLES);
#endif
			//globalLightColorResult = vec3(0.0);

		}
#endif LIGHTING


		HitResult transparentRayHit;
		dvec3 transparentRayResult = clearColor;
		double kr = 1.0;
//If it's transparent
		if (glm::epsilonNotEqual(mat.getTransparency(minRayResult.uv), 0.0, glm::epsilon<double>())) {
			bool internalOnly = false;

			dvec3 refractionRayDir = getRefractionRay(glm::normalize(minRayResult.normal), glm::normalize(ray.direction), mat.getNI(minRayResult.uv), entering, internalOnly);
			dvec3 refractionRayPos = minRayResult.position + (minRayResult.normal * (entering ? -1.0 : 1.0) * bias);

			Ray refractionRay(refractionRayPos, refractionRayDir);
			if (!internalOnly) {
				kr = shlicksApprox(mat.getNI(minRayResult.uv), minRayResult.normal, ray.inverseDirection, entering);
			}
			rayTrace(refractionRay, fi, transparentRayResult, transparentRayHit, mat.getNI(minRayResult.uv), layer + 1);
		}

		
		//TODO: this blending is not correct
		double trans = mat.getTransparency(minRayResult.uv);

#ifdef LIGHTING
		colorResult = ((1.0 - trans) * (directLightColorResult + globalLightColorResult)) + ((1.0 - (1.0 - trans)) * transparentRayResult);
#else
		colorResult = ((1.0- trans) * mat.getColor(minRayResult.uv)) + ((1.0-(1.0 - trans))*transparentRayResult);
#endif
		//colorResult = ((1.0 - trans) * dvec3(minRayResult.depth/10.0)) + ((1.0 - (1.0 - trans)) * transparentRayResult);


		return true;
	}
	colorResult = clearColor;
	return false;
}





void addModel(vector<Shape*>& shapes, string modelName, dvec3 pos = dvec3(0.0), dvec3 rot = dvec3(0.0, 0.0, 0.0)) {
	Model m(modelName, pos, rot);
	shapes.insert(shapes.end(), m.children.begin(), m.children.end());
}


AABB redoAABBs(FrameInfo& fi) {
	AABB toReturn(dvec3(0.0), dvec3(0.0));
	for (Shape* s : fi.shapes) {
		s->redoAABB(fi.currentTime);
		toReturn.encompass(s->boundingBox);
	}
	return toReturn;
}

int main()
{


	/*dvec3 q(1.0, 0.0, 0.0);//straight right
	dvec3 w(1.0, 1.0, 0.0);//diagonal up right
	dvec3 e(0.0, 1.0, 0.0);//up
	dvec3 r(-1.0, 1.0, 0.0);//diagonal up left

	printf("%f\n", glm::dot(e, q));
	printf("%f\n", glm::dot(w, q));
	printf("%f\n", glm::dot(r, q));
	printf("%f\n", glm::dot(q, q));


	exit(-1);*/


	//srand(static_cast <unsigned> (time()));

	srand(0u);

	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SRGB_CAPABLE, 1);
	glfwWindowHint(GLFW_SAMPLES, 16);
	GLFWwindow* window = glfwCreateWindow(GLsizei(frameX), GLsizei(frameY), "Renderer", NULL, NULL);
	if (window == NULL)
	{
		cout << "Window creation failed" << endl;
		exit(-1);
	}
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		cout << "GLAD init failed" << endl;
		exit(-1);
	}
	frameBufferSizeCallback(window, frameX, frameY);

	glEnable(GL_FRAMEBUFFER_SRGB);
	unsigned int VBO, VAO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glClearColor(clearColor.x, clearColor.y, clearColor.z, 1.0f);

	Shader shader("vert.glsl", "frag.glsl");
	shader.use();

	//make framebuffers

	frameBuffer = new vec3[frameX * frameY]();
	drawBuffer = new vec3[frameX * frameY]();

	//initialize textured that the framebuffer gets written to display it on a triangle
	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


	vector<Shape*> shapes;

	



	//TODO: two spheres makes it hard. NEED to do the octtree


	//Material checkers("PlainWhiteTees");
	Material checkers("PlainWhiteTees", dvec3(0.0, 0.0, 0.0), ryCheckers10x10);

	Material white("PlainWhiteTees");
	Material red("Red");
	Material green("Green");


	//Material pyt("PlainWhiteTees");


	//bottom side
	/*Vertex a(dvec3(-10.0, -10, -10.0), dvec3(0.0, 1.0, 0.0), dvec2(0.0, 0.0));
	Vertex b(dvec3(10.0, -10, -10.0), dvec3(0.0, 1.0, 0.0), dvec2(1.0, 0.0));
	Vertex c(dvec3(10.0, -10, 10.0), dvec3(0.0, 1.0, 0.0), dvec2(1.0, 1.0));
	Vertex d(dvec3(-10.0, -10, 10.0), dvec3(0.0, 1.0, 0.0), dvec2(0.0, 1.0));



	//right side
	Vertex aR(dvec3(10.0, 10, -10.0), dvec3(-1.0, 0.0, 0.0), dvec2(0.0, 1.0));
	Vertex bR(dvec3(10.0, -10, -10.0), dvec3(-1.0, 0.0, 0.0), dvec2(0.0, 0.0));
	Vertex cR(dvec3(10.0, -10, 10.0), dvec3(-1.0, 0.0, 0.0), dvec2(1.0, 0.0));
	Vertex dR(dvec3(10.0, 10, 10.0), dvec3(-1.0, 0.0, 0.0), dvec2(1.0, 1.0));

	//left side
	Vertex aL(dvec3(-10.0, 10, -10.0), dvec3(1.0, 0.0, 0.0), dvec2(1.0, 1.0));
	Vertex bL(dvec3(-10.0, -10, -10.0), dvec3(1.0, 0.0, 0.0), dvec2(1.0, 0.0));
	Vertex cL(dvec3(-10.0, -10, 10.0), dvec3(1.0, 0.0, 0.0), dvec2(0.0, 0.0));
	Vertex dL(dvec3(-10.0, 10, 10.0), dvec3(1.0, 0.0, 0.0), dvec2(0.0, 1.0));*/

	//left side



	//shapes.push_back(new Triangle(a, b, c, checkers));
	/*shapes.push_back(new Triangle(a, c, b, checkers));
	shapes.push_back(new Triangle(a, d, c, checkers));
	shapes.push_back(new Triangle(aR, bR, cR, checkers));
	shapes.push_back(new Triangle(aR, cR, dR, checkers));
	shapes.push_back(new Triangle(aL, cL, bL, checkers));
	shapes.push_back(new Triangle(aL, dL, cL, checkers));*/


	double wallRadius = 100000.0;
	double roomSize = 20.0;
	double totalRad = wallRadius + roomSize;

	shapes.push_back(new Sphere(dvec3(0.0, totalRad, 0.0), wallRadius, white, noMovement));
	shapes.push_back(new Sphere(dvec3(0.0, -totalRad, 0.0), wallRadius, white, noMovement));
	shapes.push_back(new Sphere(dvec3(totalRad, 0.0, 0.0), wallRadius, red, noMovement));
	shapes.push_back(new Sphere(dvec3(-totalRad, 0.0, 0.0), wallRadius, green, noMovement));
	//shapes.push_back(new Sphere(dvec3(0.0, 0.0, totalRad), wallRadius, Material("PlainWhiteTees"), noMovement));
	shapes.push_back(new Sphere(dvec3(0.0, 0.0, -totalRad), wallRadius, white, noMovement));

	shapes.push_back(new Sphere(dvec3(0.0, 0.0, 0.0), 4.0, Material("Mirror"), noMovement));

	shapes.push_back(new Sphere(dvec3(0.0, 69.0, 0.0), 50, Material("PlainWhiteTees", dvec3(0.5, 0.5, 0.5)), noMovement));


	//shapes.push_back(new Sphere(dvec3(1.0, -4.2, 0.0), 2.0, Material("Glass"), noMovement));

	constexpr double mypifornow = glm::pi<double>();
	//addModel(shapes, "brick2");
	//addModel(shapes, "bunny", vec3(-0.0, 0.0, -0.0));


	//shapes.push_back(new Sphere(dvec3(0.0, 0, 0), 1.5, checkers, noMovement));

	//addModel(shapes, "bunny", vec3(0.0, 0.0, 0.0), dvec3(mypifornow/2.0, 0.0, mypifornow/2.0));
	//addModel(shapes, "backpack");
	//shapes.push_back(new Sphere(dvec3(0, 0, 0), 1.0, Material("Bug"), noMovement));
	//shapes.push_back(new Sphere(dvec3(1, 1, -3), 1.0, Material("PlainWhiteTees"), noMovement));
	//shapes.push_back(new Sphere(dvec3(-1, -1, 3), 1.0, Material("PlainWhiteTees"), noMovement));

	vector<Light*> lights;

	/*lights.push_back(new PointLight(
		dvec3(5.0, 5.0, 5.0),
		dvec3(1.0, 1.0, 1.0),
		dvec3(1.0, 0.09, 0.032)
	));*/

	/*lights.push_back(new DirectionalLight(
		dvec3(-1.0, -1.0, 0.0),
		dvec3(1.0, 1.0, 1.0)
	));*/

	//lights.push_back(new SquareLight(dvec3(0.0, 10.0, 0.0), dvec3(0.0, -1.0, 0.0), 1.0, 1.0));




	FrameInfo fi;
	fi.currentTime = 0;
	fi.shapes = shapes;
	fi.lights = lights;
	fi.kdTree = NULL;


	AABB sceneBounding = redoAABBs(fi);
	printf("scene bounding min: %s, max:%s\n", glm::to_string(sceneBounding.min).c_str(), glm::to_string(sceneBounding.max).c_str());
	fi.kdTree = buildKDTree(fi.shapes, sceneBounding);

	//printKDTree(fi.kdTree);


	int frameCounter = -1;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	uint32_t fps = 12;


	time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/" + saveFileDirectory).c_str());
#ifdef OUTPUTFRAMES

	time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/"+saveFileDirectory).c_str());


	for (uint32_t frameCounter = 0; frameCounter < OUTPUTFRAMES;) {
		double currentFrame = double(frameCounter) / double(fps);
#else

	while (!glfwWindowShouldClose(window)) {
#ifdef EVERYFRAME

		double currentFrame = double(frameCounter) / double(fps);
#else
		double currentFrame = glfwGetTime();
#endif

#endif

		frameCounter++;


		fi.currentTime = currentFrame;
		deltaTime = currentFrame - lastFrame;
		//printf("that frame took %f seconds\n", deltaTime);
		lastFrame = currentFrame;

		if (int(currentFrame) > lastSecondFrameCount) {
			lastSecondFrameCount = int(currentFrame);
			float sum = 0;
			for (float f : frameTimes) {
				sum += f;
			}
			printf("fps: %f\n", sum / 30.0f);
		}
		frameTimes[frameCounter % 30] = 1.0f / float(deltaTime);
		frameTimes[frameCounter % 30] = 1.0f / float(deltaTime);

		constexpr double mypi = glm::pi<double>();
		clearBuffers();

		//dvec3 eye = vec3(sin(currentFrame/3.0) * 25.0, 0.0, cos(currentFrame/3.0) * 25.0);
		dvec3 eye = dvec3(0, 0.0, 50.0);
		dvec3 lookat = vec3(0.0, 0.0, 0.0);
		dvec3 camForward = glm::normalize(lookat - eye);
		dvec3 camUp = glm::normalize(vec3(0.0, 1.0, 0.0));
		dvec3 camRight = glm::cross(camForward, camUp);
		camUp = glm::cross(camRight, camForward);

		fi.camPosition = eye;


		double viewPortHeight = 2.0f;
		double viewPortWidth = viewPortHeight * frameRatio;

		double fov = 90;
		double focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));
		//qua rotQuat = glm::rotation(dvec3(0.0, 0.0, -1.0), camForward);




		for (int i = 0; i < MONTE_CARLO_SAMPLES; i++) {
			//clearDrawBuffer();
			if (i % 1 == 0) {
				printf("on montecarlo %i\n", i);
			}
#ifdef CONCURRENT_FOR
			concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
#else
			for (uint64_t i = 0; i < frameX * frameY; i++) {
#endif
				uint32_t x = i % frameX;
				uint32_t y = i / frameX;

		        //prd = (x == 22 && y == 103);
				double normalizedX = (double(x) / double(frameX)) - 0.5;
				double normalizedY = (double(y) / double(frameY)) - 0.5;

				dvec3 coordOnScreen = (normalizedX * camRight) + (normalizedY * camUp) + eye + (camForward * focal);
				dvec2 clipSpacePixelSize(dvec2(1.0 / double(frameX - 1.0), 1.0 / double(frameY - 1.0)));

#ifdef PIXEL_MULTISAMPLE_N
				uint32_t n = PIXEL_MULTISAMPLE_N;
#else
				uint32_t n = 1;
#endif
				dvec3 colorAcum(0);
				HitResult minRayResult;//TODO: use this for drawing the lines if you want to do that at some point
									   //TODO: change the multisampling to do a few circles instead of a square? or does it really matter, idk
				for (uint32_t x = 1; x <= n; x++) {
					double offsetX = x * (clipSpacePixelSize.x / (n + 1));
					for (uint32_t y = 1; y <= n; y++) {
						double offsetY = y * (clipSpacePixelSize.y / (n + 1));
						dvec3 rayVector = glm::normalize((coordOnScreen + dvec3(offsetX, offsetY, 0.0)) - eye);
						Ray initialRay(eye, rayVector);
						dvec3 colorOut = pathTrace(initialRay, fi);
						if(prd)printf("color out: %s\n", glm::to_string(colorOut).c_str());

						colorAcum += colorOut;


					}
				}
				drawBuffer[x + (y * frameX)] += colorAcum / double(n * n);
				if (prd) {
					drawBuffer[x + (y * frameX)] = dvec3(0.0, 0.0, 1.0);

				}
#ifdef CONCURRENT_FOR
				});
#else 
			}
#endif 
			for (uint64_t j = 0; j < frameX * frameY; j++) {
				frameBuffer[j] = drawBuffer[j] / float(i+1);
			}
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glDrawArrays(GL_TRIANGLES, 0, 3);
			glfwSwapBuffers(window);
			processInput(window);
			glfwPollEvents();
			//cin.get();
			if (i % 5 == 0) {
				saveImage((std::to_string(i/1000) + ".png"), window);
			}
		}




#ifdef OUTPUTFRAMES
		saveImage((std::to_string(frameCounter) + ".png"), window);
#endif
#ifdef CIN
		cin.get();
#endif

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}


