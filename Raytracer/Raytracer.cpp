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
/*#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image_write.h"*/



#include "Shader.hpp"
#include "Model.hpp"
//#include "WorkerPool.hpp"
#include "Triangle.hpp"
#include "Ray.hpp"


using namespace std;
using namespace glm;


uint32_t frameX = 500;
uint32_t frameY = 500;

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
Vertex transformIt(const mat4& view, const mat4& model, Vertex v) {
	//mat4 mvp = projection * view * model; //for rasterization
	mat4 mvp = view * model; //for raytracing
	mat4 normalMat = glm::inverse(glm::transpose(model));
	vec4 newPos = mvp* vec4(v.position, 1.0f);
	v.position = newPos / newPos.w;// vec3(newPos.x, newPos.y, newPos.z) / newPos.w;
	v.normal = normalMat * vec4(v.normal, 1.0f);
	return v;
}

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

void doTriangles(Mesh& mesh, const mat4& view, vector<Triangle>& triangles) {
	vector<Vertex> transformedVertices(mesh.indices.size());
	for (uint64_t i = 0; i < mesh.indices.size(); i++) {
		transformedVertices[i] = transformIt(view, mesh.getModelMat(), mesh.vertices[mesh.indices[i]]);
	}

	//construct triangles out of vertices
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
}



float angle(vec3 a, vec3 b) {
	return glm::acos(glm::dot(a, b));
}


fvec2 frameToNDC(ivec2 fc) {
	return (fvec2(floatDiv(fc.x, frameX), floatDiv(fc.y, frameY)) * 2.0f) - fvec2(1.0f);
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
	//depthBuffer = new float[frameX * frameY]();

	//initialize textured that the framebuffer gets written to display it on a triangle
	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	//Model backpack("brick");
	Model backpack = Model("backpack");

	int frameCounter = -1;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	//WorkerPool pool = WorkerPool();

	while (!glfwWindowShouldClose(window))
	{
		frameCounter++;

		//printf("starting frame %i\n", frameCounter);

		//frametime calculation
		double currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		//fps calculation
		if (int(currentFrame) > lastSecondFrameCount) {
			lastSecondFrameCount = int(currentFrame);
			float sum = 0;
			for (float f : frameTimes) {
				sum += f;
			}
			printf("fps: %f\n", sum/30.0f);
		}
		frameTimes[frameCounter % 30] = 1.0f/float(deltaTime);

		printf("currentFrame: %f\n", currentFrame);
		//projection/view matrix 
		//printf("doing matrcies\n");;
		//mat4 projection = glm::perspective(radians(70.0f), ratio, near, far);
		vec3 eye = vec3(sin(currentFrame)*5.0f, 0.0f, cos(currentFrame)*5.0f);
		//vec3 eye = vec3(0.0f, 0.0f, currentFrame);// vec3(sin(currentFrame) * 5.0f, 0.0f, cos(currentFrame) * 5.0f);

		mat4 view = glm::lookAt(eye, vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));

		//printf("clearing buffers\n");
		clearBuffers();

		vector<Triangle> triangles;

		//printf("drawing some mesh\n");
		for (Mesh& m : backpack.children) {
			doTriangles(m, view, triangles);
			//drawMesh(m, projection, view, m.transform);
		}
		printf("did triangles\n");

		float ratio = float(frameX) / float(frameY);
		float fovY = radians(70.0f);
		float fovX = fovY * ratio;
		vec3 intoScreen(0.0f, 0.0f, -1.0f);


		/*triangles = vector<Triangle>();
		float z = -50.0f;
		Vertex a = transformIt(view, mat4(1.0f), Vertex(vec3(2.0f, 0.0f, 0.0f), vec2(1.0f, 0.0f)));
		Vertex b = transformIt(view, mat4(1.0f), Vertex(vec3(0.0f, 0.0f, 2.0f), vec2(0.0f, 0.0f)));
		Vertex c = transformIt(view, mat4(1.0f), Vertex(vec3(0.0f, 2.0f, 0.0f), vec2(0.0f, 1.0f)));
		Triangle t = Triangle(a, c, b, NULL);
		triangles.push_back(t);*/


		//printf("triangles transformed, starting drawing\n");

		constexpr float mypi = glm::pi<float>();

		/*printf("triangle:\n\t%s\n\t%s\n\t%s\nbounded by:\n\t%s\n\t%s\n\n",
			glm::to_string(t[0].position).c_str(),
			glm::to_string(t[1].position).c_str(),
			glm::to_string(t[2].position).c_str(),
			glm::to_string(t.minBounding).c_str(),
			glm::to_string(t.maxBounding).c_str()
		);*/


		//printf("eye: %s\n", glm::to_string(eye).c_str());

		vec3 bary;
		vec3 boxCoord;
		vec3 hitPoint;




		for (uint64_t i = 0; i < frameX * frameY; i++) {
			{
				uint32_t x = i % frameX;
				uint32_t y = i / frameX;
		//for (uint32_t x = 0; x < frameX; x++) {
			//for (uint32_t y = 0; y < frameY; y++) {

				if (y == 0) { printf("on pixel %i, %i\n", x, y); }
				fvec2 clipCoords = frameToNDC(ivec2(x, y));
				vec3 rayVector = glm::normalize(vec3(clipCoords * mypi, -mypi));
				Ray ray(vec3(0), rayVector);
				float hitDepth = INFINITY;
				for (const Triangle& tri: triangles) {
					if (rayHit(ray, tri.minBounding, tri.maxBounding, boxCoord)) {
						if (rayHit(ray, tri, bary, hitPoint)) {
							float distance = glm::length2(boxCoord);
							if (distance < hitDepth) {
								vec2 fragUV = bary.x * tri[0].texCoords + bary.y * tri[1].texCoords + bary.z * tri[2].texCoords;
								vec3 fragColor = tri.mesh->diffuse.sample(fragUV);
								frameBuffer[x + (y * frameX)] = fragColor;
								hitDepth = distance;
							}
						}

					}
				}
			}
		}

		//printf("drew all triangles supposedly\n");
		//printf("drew all mesh\n");
		//printf("has %i triangles\n", triangles.size());


		//TODO: just put them all in one queue then have all the threads process it, instead of going as it's going

		/*queue <fragRef*> fragQueue;
		for (Triangle& tri : triangles) {
			if (tri.onscreen) {
				for (uint32_t fragx = tri.minFrags.x; fragx <= tri.maxFrags.x; fragx++) { //Lines that try to limit the box of fragments to where the triangle should be
					for (uint32_t fragy = tri.minFrags.y; fragy <= tri.maxFrags.y; fragy++) {
						fragQueue.push(new fragRef(tri, fragx, fragy));
						//pool.addJob(new fragRef(tri, fragx, fragy));
					}
				}
			}
		}*/
		//printf("gibbing que %i\n", fragQueue.size());
		//pool.gibQue(fragQueue);

		//printf("waiting for empty queue\n");
		//pool.waitForEmptyQueue();
		//printf("empty queue\n");
		//printf("moving to opengl stuff\n");
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();
		cin.get();

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}