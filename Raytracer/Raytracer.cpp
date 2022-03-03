#include <iostream>
#include <vector>
#include <fstream>
#include <mutex>
#include <execution>
#include <ppl.h>
#include <thread>
#include <filesystem>
#include <ctime>

//#include <queue>
//#include <execution>



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
//#include "Model.hpp"
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
#include "Model.hpp"

#include "Light.hpp"
#include "PointLight.hpp"
#include "DirectionalLight.hpp"


using namespace std;
using namespace std::filesystem;
using namespace glm;

#define MAXBOUNCES 6u
//#define OUTPUTFRAMES 151
#define EVERYFRAME INFINITY
#define CONCURRENT_FOR
#define KDTRACE
//#define CIN
#define LIGHTING


bool prd = false; //print debuging for refraction

uint32_t frameX = 1000;
uint32_t frameY = 1000;
double frameRatio = double(frameX) / double(frameY);

//dvec2 pixelOffset = vec2(0.5 / float(frameX), 0.5 / float(frameY));

/*float screenRatio = float(frameX) / float(frameY);
uint64_t screenY = glm::max(1000u, frameX);// scaleFactor* frameY;
uint64_t screenX = uint64_t(glm::max(1000u, frameX) * screenRatio);// scaleFactor* frameX;*/

//const float near = 0.1f;
//const float far = 50.0f;

dvec3 clearColor(0.21, 0.78, 0.95);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame
string saveFileDirectory = "";


constexpr double bias = 1e-4;


struct FrameInfo {
	vector<Shape*> shapes;
	vector<Light*> lights;
	KDNode* kdTree;
	//dmat4 view;
	dvec3 camPosition;
	double currentTime;
};


//Applies the projection and view matrix to the point
/*Vertex transformIt(const mat4& view, const mat4& model, Vertex v) {
	//mat4 mvp = projection * view * model; //for rasterization
	mat4 mvp = view * model; //for raytracing
	mat4 normalMat = glm::inverse(glm::transpose(model));
	vec4 newPos = mvp* vec4(v.position, 1.0f);
	v.position = newPos / newPos.w;// vec3(newPos.x, newPos.y, newPos.z) / newPos.w;
	v.normal = normalMat * vec4(v.normal, 1.0f);
	return v;
}*/

void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

}

vec3* frameBuffer;
//float* depthBuffer;

void frameBufferSizeCallback(GLFWwindow* window, uint64_t width, uint64_t height) {
	glViewport(0, 0, GLsizei(width), GLsizei(height));
}

void clearBuffers() {
	for (uint64_t i = 0; i < frameX * frameY; i++) {
		frameBuffer[i] = clearColor;
		//depthBuffer[i] = INFINITY;
	}
}

bool frontFacing(vec3 a, vec3 b, vec3 c) {
	mat3 m(a.x, b.x, c.x, a.y, b.y, c.y, 1.0f, 1.0f, 1.0f);
	return glm::determinant(m) > 0;
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


/*
constexpr double freeSpaceImpedance = 376.730;
//uses the fresnel equations to 
//need to return two values if you want to take light polarization into account
double reflectance(const dvec3& vectorI, const dvec3& vectorT, const dvec3& normal, double prevIOR, double newIOR) {

	//printf("calculating reflectance\n");
	double z1 = freeSpaceImpedance / prevIOR;
	double z2 = freeSpaceImpedance / newIOR;



	double cosI = glm::dot(vectorI, normal);
	double cosT = glm::dot(vectorT, -normal);

	//double costI = glm::cos(thetaI);
	//double costT = glm::cos(thetaT);

	//double cosI = glm::clamp(-1.0, 1.0, glm::dot(vectorI, vectorT));




	//double cosT = 69.0;//TODO: not real
	////etai is prior, etat is the one we're going into
	//if (cosI > 0) {
	//	swap(prevIOR, newIOR);
	//}

	////https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
	////https://www.scratchapixel.com/code.php?id=13&origin=/lessons/3d-basic-rendering/introduction-to-shading


	//double sinT = 0.0;//TODO:
	

	//printf("prevIOR: %f\n", prevIOR);
	//printf("newIOR: %f\n", newIOR);

	double z2CostI = z2 * cosI;
	double z2CostT = z2 * cosT;

	double z1CostI = z1 * cosI;
	double z1CostT = z1 * cosT;

	//printf("z2CostI: %f\n", z2CostI);
	//printf("z2CostT: %f\n", z2CostT);
	//printf("z1CostI: %f\n", z1CostI);
	//printf("z1CostT: %f\n", z1CostT);

	double rS = glm::pow(glm::abs((z2CostI - z1CostT) / (z2CostI + z1CostT)), 2.0);
	double rP = glm::pow(glm::abs((z2CostT - z1CostI) / (z2CostT + z1CostI)), 2.0);




	return (rS + rP) / 2.0;
}*/


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

	//TODO: navigate octtree instead
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






bool rayTrace(const Ray& ray, const FrameInfo& fi, dvec3& colorResult, HitResult& minRayResult, double currentIOR = 1.0, uint32_t layer = 0)
{



	if (layer > MAXBOUNCES) {
		colorResult = clearColor;
		return false;
	}

	colorResult = dvec3(0);
	minRayResult = shootRay(ray, fi); //ray is in camera space
	if (minRayResult.shape != NULL) {
		

		
		Material mat = minRayResult.shape->mat;

		double hitAngle = glm::acos(glm::dot(minRayResult.normal, ray.inverseDirection));

		bool entering = hitAngle < (glm::pi<double>() / 2.0);




		//TODO: apparently i messedup something? http://cosinekitty.com/raytrace/chapter09_refraction.html


		dvec3 lightColorResult = dvec3(0.0);
		//double kr = glm::clamp(shlicksApprox(currentIOR, mat.ni, minRayResult.normal, ray.inverseDirection), 0.0, 1.0);

		//printf("kr: %f\n", kr);

		for (Light* light : fi.lights) {
			Ray lightRay = light->getRay(minRayResult.position+(minRayResult.normal*bias));//should be in camera space

			HitResult lightRayResult = shootRay(lightRay, fi);
			double lightDistance = light->getDistance(minRayResult.position);
			if (lightDistance < lightRayResult.depth) {//shawows
				double attenuation = light->getAttenuation(lightDistance);
				dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal)- lightRay.direction);
				dvec3 H = glm::normalize(lightRay.direction + ray.inverseDirection);
				//double spec = glm::pow(glm::max(0.0, glm::dot(minRayResult.normal, H)), minRayResult.shape->mat.ns);//to the specular exponent
				double spec = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, ray.inverseDirection)), mat.getNS(minRayResult.uv));//to the specular exponent
				dvec3 specular = (light->color * spec) / attenuation;
				double diff = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
				dvec3 diffuse = (mat.getColor(minRayResult.uv) * light->color * diff) / attenuation;
				//dvec3 diffuse = (dvec3(0.0, 0.0, 1.0) * light->color * diff) / attenuation;
				lightColorResult += (diffuse + specular);
			}

		}


		HitResult transparentRayHit;
		dvec3 transparentRayResult = clearColor;
		double kr = 1.0;

		//bool transparentRayDidhit = false;

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
		colorResult = ((1.0 - trans) * lightColorResult) + ((1.0 - (1.0 - trans)) * transparentRayResult);
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


//TODO: implement a global sphere for a global color
int main()
{
	/*dvec3 n = dvec3(0.0, -1.0, 0.0);
	dvec3 x = dvec3(1.0, 1.0, 0.0);
	dvec3 y = dvec3(1.0, -0.1, 0.0);

	printf("%f\n", glm::dot(n, x));
	printf("%f\n", glm::dot(n, y));

	exit(-1);*/




	/*dvec3 normal = dvec3(0.0, 1.0, 0.0);
	dvec3 incident = glm::normalize(dvec3(2.0, -1.0, 0.0));




	dvec3 other = getRefractionRay(normal, incident, 1.0, 1.54);

	double inAngle = glm::acos(glm::dot(-normal, incident));
	double newangle = glm::acos(glm::dot(normal, other));

	printf("inAngle: %f\n", degrees(inAngle));
	printf("newangle: %f\n", degrees(newangle));

	exit(-1);*/
	

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


	Material checkers("PlainWhiteTees", ryCheckers10x10);



	Vertex a(dvec3(-5.0, -5, -10.0), dvec2(0.0, 0.0));
	Vertex b(dvec3(5.0, -5, -10.0), dvec2(1.0, 0.0));
	Vertex c(dvec3(5.0, -5, 10.0), dvec2(1.0, 1.0));
	Vertex d(dvec3(-5.0, -5, 10.0), dvec2(0.0, 1.0));

	//shapes.push_back(new Triangle(a, b, c, checkers));

	//shapes.push_back(new Triangle(a, c, b, checkers));
	//shapes.push_back(new Triangle(a, d, c, checkers));

	constexpr double mypifornow = glm::pi<double>();
	//addModel(shapes, "brick2");
	addModel(shapes, "bunny", vec3(-0.0, 0.0, -0.0));
	//addModel(shapes, "bunny", vec3(-5.0, 0.0, 5.0));
	//addModel(shapes, "bunny", vec3(5.0, 0.0, -5.0));
	//addModel(shapes, "bunny", vec3(5.0, 0.0, 5.0));
	//addModel(shapes, "bunny", vec3(0.0, 0.0, -5.0));
	//addModel(shapes, "bunny", vec3(0.0, 0.0, 5.0));
	//addModel(shapes, "bunny", vec3(-5.0, 0.0, 0.0));
	//addModel(shapes, "bunny", vec3(5.0, 0.0, 0.0));

	shapes.push_back(new Sphere(dvec3(5.0, 1, 0), 1.5, Material("Glass"), noMovement));

	//addModel(shapes, "bunny", vec3(0.0, 0.0, 0.0), dvec3(mypifornow/2.0, 0.0, mypifornow/2.0));
	//addModel(shapes, "backpack");
	//shapes.push_back(new Sphere(dvec3(0, 0, 0), 1.0, Material("Bug"), noMovement));
	//shapes.push_back(new Sphere(dvec3(1, 1, -3), 1.0, Material("PlainWhiteTees"), noMovement));
	//shapes.push_back(new Sphere(dvec3(-1, -1, 3), 1.0, Material("PlainWhiteTees"), noMovement));

	vector<Light*> lights;

	/*lights.push_back(new PointLight(
		dvec3(0.0, 0, 5.0),
		dvec3(1.0, 1.0, 1.0),
		dvec3(1.0, 0.09, 0.032)
	));*/

	lights.push_back(new DirectionalLight(
		dvec3(-1.0, -1.0, -1.0),
		dvec3(1.0, 1.0, 1.0)
	));




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
#ifdef EVERYFRAME

	for (uint32_t frameCounter = 0; frameCounter < INFINITY;) {
		double currentFrame = double(frameCounter) / double(fps);

#else
	while (!glfwWindowShouldClose(window)) {
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
			//printf("fps: %f\n", sum / 30.0f);
		}
		frameTimes[frameCounter % 30] = 1.0f / float(deltaTime);
		frameTimes[frameCounter % 30] = 1.0f / float(deltaTime);

		constexpr double mypi = glm::pi<double>();
		clearBuffers();

		dvec3 eye = vec3(sin(currentFrame) * 10.0, 2.0, cos(currentFrame) * 10.0);
		//dvec3 eye = dvec3(0, 2.0, 10.0);
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
		qua rotQuat = glm::rotation(dvec3(0.0, 0.0, -1.0), camForward);





#ifdef CONCURRENT_FOR
		concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
#else
		for (uint64_t i = 0; i < frameX * frameY; i++) {
#endif
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;

			//prd =(x == 250 && y == 140);
			


			double normalizedX = (double(x) / double(frameX)) -0.5;
			double normalizedY = (double(y) / double(frameY)) -0.5;

			dvec3 coordOnScreen = (normalizedX * camRight) + (normalizedY * camUp) + eye + (camForward*focal);
			dvec2 clipSpacePixelSize(dvec2(1.0 / double(frameX - 1.0), 1.0 / double(frameY - 1.0)));

			//multisample n*n
			uint32_t n = 1;
			dvec3 colorAcum(0);

			HitResult minRayResult;//TODO: use this for drawing the lines if you want to do that at some point
								   //TODO: change the multisampling to do a few circles instead of a square? or does it really matter, idk

			for (uint32_t x = 1; x <= n; x++) {
				double offsetX = x * (clipSpacePixelSize.x / (n + 1));
				for (uint32_t y = 1; y <= n; y++) {
					double offsetY = y * (clipSpacePixelSize.y / (n + 1));

					dvec3 rayVector = glm::normalize((coordOnScreen + dvec3(offsetX, offsetY, 0.0))-eye);
					dvec3 colorOut;

					Ray initialRay(eye, rayVector);
					bool hit = rayTrace(initialRay, fi, colorOut, minRayResult, 1.0, 0);
					if (hit) {
						colorAcum += colorOut;
					}
					else {
						colorAcum += clearColor;
					}

				}
			}
			frameBuffer[x + (y * frameX)] = colorAcum / double(n * n);

			if (prd) {
				frameBuffer[x + (y * frameX)] = dvec3(0.0, 1.0, 1.0);
					
			}
#ifdef CONCURRENT_FOR
		});
#else 
		}
#endif 
		//printf("did parallelForloop\n");


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();
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


