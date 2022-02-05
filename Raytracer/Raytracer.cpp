#include <iostream>
#include <vector>
#include <fstream>
#include <mutex>
#include <execution>
#include <ppl.h>
#include <thread>
//#include <queue>
//#include <execution>



//#define GLM_FORCE_LEFT_HANDED 
//#define GLM_RIGHT_HANDED 

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

#include "Light.hpp"
#include "PointLight.hpp"

using namespace std;
using namespace glm;


uint32_t framesToRender = 1;

uint32_t frameX = 250;
uint32_t frameY = 250;

fvec2 pixelOffset = vec2(0.5 / float(frameX), 0.5 / float(frameY));

float screenRatio = float(frameX) / float(frameY);
uint64_t screenY = 1000;// scaleFactor* frameY;
uint64_t screenX = uint64_t(1000.0f* screenRatio);// scaleFactor* frameX;

const float near = 0.1f;
const float far = 50.0f;

vec3 clearColor(0.0f, 0.0f, 0.0f);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame


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

//mutex perFrameMutex;
//condition_variable perFrameCondition;

/*void drawMesh(Mesh& mesh, const mat4& projection, const mat4& view, const mat4& parentTx) {
	vector<Vertex> transformedVertices(mesh.indices.size());
	for (uint64_t i = 0; i < mesh.indices.size(); i++) {
		transformedVertices[i] = transformIt(view, mesh.getModelMat(), mesh.vertices[mesh.indices[i]]);
	}


	vector<Triangle> triangles;
	//construct  out of vertices
	for (int i = 0; i < transformedVertices.size(); i += 3) {
		if (frontFacing(transformedVertices[i + 0].position,
			transformedVertices[i + 1].position,
			transformedVertices[i + 2].position)) {

			triangles.push_back(Triangle(
				transformedVertices[i + 0],
				transformedVertices[i + 1],
				transformedVertices[i + 2], 
				&mesh
			));
		}
	}


	
	for (Triangle& tri : triangles) {
		if (tri.onscreen) {
			for (uint32_t fragx = tri.minFrags.x; fragx <= tri.maxFrags.x; fragx++) { //Lines that try to limit the box of fragments to where the triangle should be
				for (uint32_t fragy = tri.minFrags.y; fragy <= tri.maxFrags.y; fragy++) {
					fvec2 fragCoord = fvec2((floatDiv(fragx, frameX) * 2.0f) - 1.0f, (floatDiv(fragy, frameY) * 2.0f) - 1.0f) + pixelOffset;
					if (pointInTriangle(fragCoord, tri)) {
						vec3 baryCoords = tri.getBaryCoords(fragCoord);
						if (baryCoords.x >= 0 && baryCoords.y >= 0 && baryCoords.z >= 0) {
							float fragDepth = glm::dot(vec3(tri[0].position.z, tri[1].position.z, tri[2].position.z), baryCoords);

							vec2 fragUV = mat3x2(tri[0].texCoords, tri[1].texCoords, tri[2].texCoords) * baryCoords;
							vec3 fragColor = tri.mesh->diffuse.sample(fragUV);
							if (fragDepth < depthBuffer[fragx + (fragy * frameX)]) {
								depthBuffer[fragx + (fragy * frameX)] = fragDepth;
								frameBuffer[fragx + (fragy * frameX)] = fragColor;
							}
						}
					}
				}
			}
		}
	}
}*/

/*void doTriangles(Mesh& mesh, const mat4& view, vector<Triangle>* triangles, vec3& minSceneBounding, vec3& maxSceneBounding) {
	vector<Vertex> transformedVertices(mesh.indices.size());
	for (uint64_t i = 0; i < mesh.indices.size(); i++) {
		transformedVertices[i] = transformIt(view, mesh.getModelMat(), mesh.vertices[mesh.indices[i]]);
	}

	//construct triangles out of vertices
	for (int i = 0; i < transformedVertices.size(); i += 3) {
		for (int k = 0; k < 3; k++) {
			minSceneBounding = glm::min(minSceneBounding, transformedVertices[i + k].position);
			maxSceneBounding = glm::max(maxSceneBounding, transformedVertices[i + k].position);
		}

		if (frontFacing(transformedVertices[i + 0].position,
			transformedVertices[i + 1].position,
			transformedVertices[i + 2].position)) {

			triangles->push_back(Triangle(
				transformedVertices[i + 0],
				transformedVertices[i + 1],
				transformedVertices[i + 2],
				&mesh
			));
		}
	}
}*/



float angle(vec3 a, vec3 b) {
	return glm::acos(glm::dot(a, b));
}


float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}

//i misunderstood what ndc was when I made this, named incorrectly
fvec2 frameToNDC(ivec2 fc) {
	return (fvec2(floatDiv(fc.x, frameX), floatDiv(fc.y, frameY)) * 2.0f) - fvec2(1.0f);
}


//since it's not part of the assignment to write a function to save to image, i just copied this directly from
//https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/
void saveImage(char* filepath, GLFWwindow* w) {
	int width, height;
	//glfwGetFramebufferSize(w, &width, &height);

	width = frameX;
	height = frameY;
	GLsizei nrChannels = 3;
	GLsizei stride = nrChannels * width;
	stride += (stride % 4) ? (4 - stride % 4) : 0;
	GLsizei bufferSize = stride * height;
	std::vector<char> buffer(bufferSize);
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, frameBuffer);
	stbi_flip_vertically_on_write(true);
	stbi_write_png(filepath, frameX, frameY, nrChannels, frameBuffer, stride);
}

void printVector(vector<int> v) {
	printf("[");
	for (const int& x : v) {
		printf("%i, ", x);
	}
	printf("]\n");
}





HitResult shootRay(const Ray& ray, const vector<Shape*>& shapes, const dmat4& view, double currentTime) {

	HitResult minRayResult;
	minRayResult.shape = NULL;
	minRayResult.depth = INFINITY;
	for (Shape* shape : shapes) {
		HitResult rayResult;
		//if (shape->rayAABB(ray, view)) {
			if (shape->rayHit(ray, rayResult, view, currentTime)) {
				if (rayResult.depth < minRayResult.depth) {
					minRayResult = rayResult;
				}
			}
		//}
	}
	return minRayResult;
}

bool rayTrace(const Ray& ray, const vector<Shape*>& shapes, const vector<Light*>& lights, const mat4& view, vec3& colorResult, double currentTime)
{
	double bias = 1e-4;
	//constexpr float epsilon = glm::epsilon<float>();
	colorResult = dvec3(0);
	HitResult minRayResult = shootRay(ray, shapes, view, currentTime);
	if (minRayResult.shape != NULL) {
		//printf("hit position: %s\n", glm::to_string(minRayResult.position).c_str());
		dvec3 diffuse = dvec3(0.0);

		for (Light* light : lights) {

			Ray lightRay = light->getRay(minRayResult.position+(minRayResult.normal*bias), view);
			HitResult lightRayResult = shootRay(lightRay, shapes, view, currentTime);

			double lightDistance = light->getDistance(minRayResult.position, view);

			//printf("light distance2 %f\n", lightDistance2);
			if (lightDistance < lightRayResult.depth) {

				double diff = glm::max(glm::dot(minRayResult.normal, lightRay.direction), 0.0);
				//printf("supposed distance to light: %f\n", distanceToLight);
				//printf("actual distance to light: %f\n", glm::sqrt(distanceToLight));
				diffuse += light->color * diff / light->getAttenuation(lightDistance);

				//diffuse = dvec3(lightDistance2 / 50.0);
				//printf("distance to light: %f\n", sqrt(distanceToLight));

				//take the light into account, otherwise don't add the light value
			}
		}


		colorResult = diffuse * minRayResult.shape->getColor(minRayResult);
		//colorResult = minRayResult.shape->getColor(minRayResult);


		return true;
	}
	return false;


}

int main()
{
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SRGB_CAPABLE, 1);
	glfwWindowHint(GLFW_SAMPLES, 16);
	GLFWwindow* window = glfwCreateWindow(GLsizei(screenX), GLsizei(screenY), "Renderer", NULL, NULL);
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
	frameBufferSizeCallback(window, screenX, screenY);

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
	shapes.push_back(new Sphere(dvec3(0), 1.0, oscilateX));
	shapes.push_back(new Sphere(dvec3(0, -5, 0), 3.0, noMovement));
	shapes.push_back(new Sphere(dvec3(3, 0, 0), 1.0, noMovement));

	vector<Light*> lights;

	lights.push_back(new PointLight(
		dvec3(1.0, 5.0, 0.0),
		dvec3(1.0, 1.0, 1.0),
		dvec3(1.0, 0.09, 0.032)
	));




	uint32_t fps = 24;
	while (!glfwWindowShouldClose(window)) {
	//for (uint32_t frameCounter = 0; frameCounter < framesToRender; frameCounter++) {
		printf("working frame %i\n", frameCounter);
		//float currentFrame = float(frameCounter) / float(fps);

		double currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		//((PointLight*)lights.at(0))->position = vec3(sin(currentFrame * 2.0f) * 20.0f, 5.0f, 0.0f);

		constexpr float mypi = glm::pi<float>();
		//vec3 eye = vec3(sin(currentFrame) * 10.0f, 2.0f, cos(currentFrame) * 10.0f);
		vec3 eye = vec3(0.0f, 2.0f, 5.0f);

		mat4 view = glm::lookAt(eye, vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));

		//printf("clearing buffers\n");
		clearBuffers();

		//vector<Triangle>* triangles = new vector<Triangle>();

		//for (uint64_t i = 0; i < frameX * frameY; i++) {
		concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;

			vec2 clipCoords = frameToNDC(ivec2(x, y));
			vec3 rayVector = glm::normalize(vec3(clipCoords * mypi, -mypi));//TODO: this is incorrect
			Ray initialRay(vec3(0), rayVector);
			vec3 colorResult;

			bool hit = rayTrace(initialRay, shapes, lights, view, colorResult, currentFrame);
			if (hit) {
				frameBuffer[x + (y * frameX)] = colorResult;
			}
		});
		//}

		printf("did parallelForloop\n");


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();

		saveImage((char*)("out/"+std::to_string(frameCounter) + ".png").c_str(), window);
		//cin.get();

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}


