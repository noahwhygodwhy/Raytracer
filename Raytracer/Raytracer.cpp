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
#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image_write.h"




#include "Shader.hpp"
//#include "Model.hpp"
//#include "WorkerPool.hpp"
#include "Triangle.hpp"
#include "Ray.hpp"
#include "KDNode.hpp"
#include "Shape.hpp"
#include "Sphere.hpp"
#include "Biconvex.hpp"

#include "Light.hpp"
#include "PointLight.hpp"
#include "DirectionalLight.hpp"


using namespace std;
using namespace std::filesystem;
using namespace glm;

#define MAXBOUNCES 8
//#define OUTPUTFRAMES 151
#define CONCURRENT_FOR



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


struct FrameInfo {
	vector<Shape*> shapes;
	vector<Light*> lights;
	dmat4 view;
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



//uses the fresnel equations to 
//need to return two values if you want to take light polarization into account
/*double reflectance(const dvec3& vectorI, const dvec3& vectorT, double z1, double z2, double prevIOR, double newIOR) {





	//double costI = glm::cos(thetaI);
	//double costT = glm::cos(thetaT);

	double cosI = glm::clamp(-1.0, 1.0, glm::dot(vectorI, vectorT));

	double cosT = 69.0;//TODO: not real
	//etai is prior, etat is the one we're going into
	if (cosI > 0) {
		swap(prevIOR, newIOR);
	}

	//https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
	//https://www.scratchapixel.com/code.php?id=13&origin=/lessons/3d-basic-rendering/introduction-to-shading


	double sinT = 0.0;//TODO:

	double z2CostI = z2 * cosI;
	double z2CostT = z2 * cosT;

	double z1CostI = z1 * cosI;
	double z1CostT = z1 * cosT;

	double rS = glm::pow(glm::abs((z2CostI - z1CostT) / (z2CostI - z1CostT)), 2.0);
	double rP = glm::pow(glm::abs((z2CostT - z1CostI) / (z2CostT - z1CostI)), 2.0);

	return (rS + rP) / 2.0;
}*/


double shlicksApprox(double prevIOR, double newIOR, const vec3& normal, const vec3& viewVector) {
	double r0 = glm::pow((prevIOR - newIOR) / (prevIOR + newIOR), 2.0);

	double rTheta = r0 + ((1 - r0) * (1 - glm::acos(glm::dot(normal, viewVector))));
	return rTheta;
}

HitResult shootRay(const Ray& ray, const FrameInfo& fi) {

	HitResult minRayResult;
	minRayResult.shape = NULL;
	minRayResult.depth = INFINITY;

	//TODO: navigate octtree instead
	for (Shape* shape : fi.shapes) {
		HitResult rayResult;
		//if (shape->rayAABB(ray, view)) { //TODO:
			if (shape->rayHit(ray, rayResult, fi.view, fi.currentTime)) {
				if (rayResult.depth < minRayResult.depth) {
					//printf("hit something\n");
					minRayResult = rayResult;
				}
			}
		//}
	}
	return minRayResult;
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
		return false;
	}

	
	double bias = 1e-4;
	colorResult = dvec3(0);
	minRayResult = shootRay(ray, fi); //ray is in camera space
	if (minRayResult.shape != NULL) {

		//yeah honestly fuck this, do this instead? https://www.scratchapixel.com/lessons/3d-basic-rendering/global-illumination-path-tracing/global-illumination-path-tracing-practical-implementation
		
		
		/*if (glm::epsilonNotEqual(minRayResult.shape->mat.opaqueness, 1.0, glm::epsilon<double>())) {
			//TODO: do stuff for transparent reflective object
		}
		else if (glm::epsilonNotEqual(minRayResult.shape->mat.metalness, 0.0, glm::epsilon<double>())) {

			double roughness = 1.0 - minRayResult.shape->mat.metalness;


			dvec3 reflectionColor = castRay(hitPoint, reflectionDirection, ..., depth + 1);
			dvec3 refractionColor = castRay(hitPoint, refractionDirection, ..., depth + 1);
			

			double reflection = reflectance(ray.direction, minRayResult.normal, currentIOR, minRayResult.shape->mat.ni);
			colorResult = (reflectionColor * reflection) + (refractionColor * (1.0 - reflection));


			dvec3 colorFromReclection = dvec3(0.0);
			dvec3 colorFromAlbedo = minRayResult.shape->mat.kd;



			//TODO: do stuff for reflective object
		}
		else {
			//TODO: do stuff for phong object
		}*/



		//dvec3 diffuse = dvec3(0.0);
		//dvec3 specular = dvec3(0.0);
		//dvec3 ambient = (clearColor * 0.1)* minRayResult.shape->mat.ka;//TODO this is wrong, like really, be bettter

		Material mat = minRayResult.shape->mat;

		double hitAngle = glm::acos(glm::dot(minRayResult.normal, ray.inverseDirection));

		bool entering = hitAngle < (glm::pi<double>() / 2.0);


		HitResult transparentRayHit;
		dvec3 transparentRayResult = clearColor;
		if (mat.transparency > 0.0) {
			//based on https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
			dvec3 ncs = glm::cross(minRayResult.normal, ray.inverseDirection);
			double inSqrt = 1.0-((glm::pow(currentIOR / mat.ni, 2.0) * glm::dot(ncs, ncs)));
			dvec3 rightSide = minRayResult.normal * glm::sqrt(inSqrt);
			dvec3 leftSide = (currentIOR / mat.ni) * glm::cross(minRayResult.normal, glm::cross(-minRayResult.normal, ray.inverseDirection));

			dvec3 refractionRayDir = glm::normalize(leftSide - rightSide);
			dvec3 refractionRayPos = minRayResult.position + (minRayResult.normal * (entering ? -1.0 : 1.0));

			Ray refractionRay(refractionRayPos, refractionRayDir);

			rayTrace(ray, fi, transparentRayResult, transparentRayHit, mat.ni, layer + 1);

		}




		//colorResult += dvec3(0.05)*mat.kd;

		dvec3 lightColorResult = dvec3(0.0);
		double kr = shlicksApprox(currentIOR, mat.ni, minRayResult.normal, ray.inverseDirection);

		//printf("kr: %f\n", kr);

		for (Light* light : fi.lights) {

			Ray lightRay = light->getRay(minRayResult.position+(minRayResult.normal*bias), fi.view);//should be in camera space

			//printf("light ray origin: %s, direction: %s\n", glm::to_string(lightRay.origin).c_str(), glm::to_string(lightRay.direction).c_str());

			HitResult lightRayResult = shootRay(lightRay, fi);


			double lightDistance = light->getDistance(minRayResult.position, fi.view);
			if (lightDistance < lightRayResult.depth) {//shawows



				double attenuation = light->getAttenuation(lightDistance);

				//bad attempt at porting strauss to the raytracer
				//dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal) - lightRay.direction);
				//double specular = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, ray.inverseDirection)), minRayResult.shape->mat.ns);//to the specular exponent
				//double diffuse = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
				//double attenuation = light->getAttenuation(lightDistance);
				//double angleOfIncidence = cos(dot(minRayResult.normal, lightReflectVector));
				//dvec3 vecToCamera = fi.camPosition - minRayResult.position;
				//double viewAngle = cos(dot(minRayResult.normal, vecToCamera));
				//double diffuseReflectivity = (1.0 - (pow(mat.smoothness, 3.0))) * (1.0 - mat.transparency);
				//double diffuseReflectivityMultiplier = (1.0 - (mat.metalness * mat.smoothness));
				//double rn = (1.0 - mat.transparency) - diffuseReflectivity;
				//double j = F(angleOfIncidence) * G(angleOfIncidence) * G(viewAngle);
				//double adjsutedReflectivity = glm::min(1.0, rn + ((rn + kj) * j));
				//double specularReflectivity = specular * adjsutedReflectivity;
				//dvec3 specularColor = dvec3(1.0) + (mat.metalness * (1.0 - F(angleOfIncidence)) * (mat.kd - dvec3(1.0)));
				//dvec3 specularContrib = specularReflectivity * specularColor;
				//dvec3 diffuseContrib = diffuse * diffuseReflectivityMultiplier * diffuseReflectivity * mat.kd;
				//colorResult += (diffuseContrib + specularContrib) * attenuation;




				
				dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal)- lightRay.direction);

				//ambient += minRayResult.shape->mat.ka;


				dvec3 H = glm::normalize(lightRay.direction + ray.inverseDirection);
				double spec = glm::pow(glm::max(0.0, glm::dot(minRayResult.normal, H)), minRayResult.shape->mat.ns);//to the specular exponent

				//double spec = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, ray.inverseDirection)), minRayResult.shape->mat.ns);//to the specular exponent
				dvec3 specular = (light->color * spec) / attenuation;
				double diff = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
				dvec3 diffuse = (minRayResult.shape->mat.kd*light->color * diff) / attenuation;

				//colorResult += Kr * reflectionColor + (1 - Kr) * refractionColor;
				lightColorResult += (diffuse + specular);

			}

		}
		
		colorResult =  (kr*lightColorResult) + ((1.0 - kr) * mat.transparency * transparentRayResult);
		

		//double kr = reflectance();//
		// Schlick’s Approximation instead of fresnel's equation?x`
		//fresnel(dir, hitNormal, isect.hitObject->ior, kr);
		//bool oddSecond = uint32_t(floor(fi.currentTime)) % 2 == 0;

		//colorResult = (ambient + diffuse + specular) * minRayResult.shape->getColor(minRayResult);
		//colorResult = minRayResult.shape->getColor(minRayResult);
		//colorResult = dvec3(minRayResult.uv, 0.0);
		//colorResult = minRayResult.shape->getColor(minRayResult);
		return true;
	}
	return false;
}



//TODO: implement a global sphere for a global color
int main()
{

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

	//Model backpack("brick");
	//Model backpack = Model("backpack");

	int frameCounter = -1;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	vector<Shape*> shapes;

	



	//TODO: two spheres makes it hard. NEED to do the octtree


	//the two whitted spheres
	//shapes.push_back(new Sphere(dvec3(2, 2, -2), 2, materials.at("Polished Silver"), noMovement));
	//shapes.push_back(new Sphere(dvec3(-2, 4, 2), 2.5, materials.at("Polished Silver"), noMovement));

	
	//four axis spheres
	//shapes.push_back(new Sphere(dvec3(0, 0, 0), 0.5, materials.at("Polished Silver"), noMovement));
	//shapes.push_back(new Sphere(dvec3(2, 0, 0), 0.5, materials.at("Polished Silver"), noMovement));
	//shapes.push_back(new Sphere(dvec3(0, 2, 0), 0.5, materials.at("Polished Silver"), noMovement));
	//shapes.push_back(new Sphere(dvec3(0, 0, 2), 0.5, materials.at("Polished Silver"), noMovement));
	

	shapes.push_back(new Sphere(dvec3(0, 2, 0), 5.0, materials.at("Glass"), noMovement));

	//shapes.push_back(new Biconvex(dvec3(0, 1, 0), dvec3(0.0, 0.0, 1.0), 2.0, 20.0, noMovement));

	/*Vertex a(dvec3(-7.0, -3.0, -30.0), dvec2(0.0, 0.0));
	Vertex b(dvec3(13.0, -3.0, -30.0), dvec2(1.0, 0.0));
	Vertex c(dvec3(13.0, -3.0, 10.0), dvec2(1.0, 1.0));
	Vertex d(dvec3(-7.0, -3.0, 10.0), dvec2(0.0, 1.0));*/

	Vertex a(dvec3(-5.0, 0, -10.0), dvec2(0.0, 0.0));
	Vertex b(dvec3(5.0, 0, -10.0), dvec2(1.0, 0.0));
	Vertex c(dvec3(5.0, 0, 10.0), dvec2(1.0, 1.0));
	Vertex d(dvec3(-5.0, 0, 10.0), dvec2(0.0, 1.0));

	//shapes.push_back(new Triangle(a, c, b, NULL));
	//shapes.push_back(new Triangle(a, d, c, NULL));

	//shapes.push_back(new Sphere(dvec3(0, 0, 0), 1.0, materials.at("Polished Silver"), noMovement));

	vector<Light*> lights;

	/*lights.push_back(new PointLight(
		dvec3(0.0, 0, 5.0),
		dvec3(1.0, 1.0, 1.0),
		dvec3(1.0, 0.09, 0.032)
	));*/

	lights.push_back(new DirectionalLight(
		dvec3(0.0, -1.0, -1.0),
		dvec3(1.0, 1.0, 1.0)
	));




	FrameInfo fi;
	fi.currentTime = 0;
	fi.shapes = shapes;
	fi.lights = lights;
	fi.view = mat4(1.0);



	uint32_t fps = 24;
#ifdef OUTPUTFRAMES

	time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/"+saveFileDirectory).c_str());

	for (uint32_t frameCounter = 0; frameCounter < OUTPUTFRAMES; frameCounter++) {
		double currentFrame = double(frameCounter) / double(fps);
#else
	while (!glfwWindowShouldClose(window)) {
		double currentFrame = glfwGetTime();
#endif
		//printf("working frame %i\n", frameCounter);

		fi.currentTime = currentFrame;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		//((PointLight*)lights.at(0))->position = vec3(sin(currentFrame * 2.0f) * 20.0f, 5.0f, 0.0f);

		constexpr double mypi = glm::pi<double>();
		//dvec3 eye = vec3(-5.0, 0.0, 10);
		clearBuffers();


		dvec3 cameraForward = dvec3(sin(currentFrame), 0.0, cos(currentFrame));
		//dvec3 cameraForward = dvec3(0.0, 0.0, -1.0);


		dvec3 cameraUp = dvec3(0.0, 1.0, 0.0);

		dvec3 eye = vec3(sin(currentFrame) * 10.0, 2.0, cos(currentFrame) * 10.0);
		//dvec3 eye = vec3(0.0, 2.0, 10.0); 
		//dmat4 view = glm::lookAt(eye, eye+cameraForward, cameraUp);
		fi.view = glm::lookAt(eye, dvec3(0.0, 3.0, 0.0), cameraUp);


		fi.camPosition = transformPos(eye, dmat4(1.0), fi.view);


		//printf("===========\n");
		//printf("reveye: %s\n", glm::to_string(glm::normalize(-eye)).c_str());
		//printf("light reversedir: %s\n", glm::to_string(((DirectionalLight*)lights.at(0))->getRay(vec3(0), fi.view).direction).c_str());



		double viewPortHeight = 2.0f;
		double viewPortWidth = viewPortHeight * frameRatio;

		double fov = 70;
		//printf("current frame: %f\n", fov);
		double focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));
		dvec3 horizontal = dvec3(viewPortWidth, 0, 0);
		dvec3 vertical = dvec3(0, viewPortHeight, 0);

		dvec2 hv(viewPortWidth, viewPortHeight);
		dvec3 origin(0);
		dvec3 lowerLeftCorner = origin - (horizontal / 2.0) - (vertical / 2.0) - dvec3(0, 0, focal);

#ifdef CONCURRENT_FOR
		concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
#else
		for (uint64_t i = 0; i < frameX * frameY; i++) {
#endif
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;

			dvec2 uv(double(x) / double((frameX - 1.0)), double(y) / double((frameY - 1.0)));
			dvec2 clipSpacePixelSize(dvec2(1.0 / double(frameX - 1.0), 1.0 / double(frameY - 1.0)));

			//multisample n*n
			uint32_t n = 1;
			dvec3 colorAcum(0);

			dvec2 pixelBottomLeft = ((uv.x * horizontal) + (uv.y * vertical)).xy - (clipSpacePixelSize / 2.0);

			HitResult minRayResult;//TODO: use this for drawing the lines if you want to do that at some point
								   //TODO: change the multisampling to do a few circles instead of a square? or does it really matter, idk

			for (uint32_t x = 1; x <= n; x++) {
				double offsetX = x * (clipSpacePixelSize.x / (n + 1));
				for (uint32_t y = 1; y <= n; y++) {
					double offsetY = y * (clipSpacePixelSize.y / (n + 1));

					dvec3 rayVector = lowerLeftCorner + (uv.x * horizontal) + (uv.y * vertical) + dvec3(offsetX, offsetY, 0.0);
					dvec3 colorOut;

					Ray initialRay(dvec3(0.0, 0.0, focal), glm::normalize(rayVector));
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
		//cin.get();

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}


