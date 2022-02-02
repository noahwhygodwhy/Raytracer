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
#include "Model.hpp"
//#include "WorkerPool.hpp"
#include "Triangle.hpp"
#include "Ray.hpp"
#include "KDNode.hpp"

using namespace std;
using namespace glm;


uint32_t framesToRender = 24;

uint32_t frameX = 200;
uint32_t frameY = 200;

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

//i misunderstood what ndc was when I made this, named incorrectly
fvec2 frameToNDC(ivec2 fc) {
	return (fvec2(floatDiv(fc.x, frameX), floatDiv(fc.y, frameY)) * 2.0f) - fvec2(1.0f);
}


//since it's not part of the assignment to write a function to save to image, i just copied this directly from
//https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/
void saveImage(char* filepath, GLFWwindow* w) {
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
	stbi_write_png(filepath, width, height, nrChannels, buffer.data(), stride);
}

void printVector(vector<int> v) {
	printf("[");
	for (const int& x : v) {
		printf("%i, ", x);
	}
	printf("]\n");
}

int main()
{
	/*vector<int> a = {0, 1, 2, 3, 4, 5, 6, 7, 8};

	vector<int> b = vector<int>(a.begin(), a.begin()+1);

	printVector(a);
	printVector(b);

	

	exit(0);*/
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

	uint32_t fps = 24;

	for (uint32_t frameCounter = 0; frameCounter < framesToRender; frameCounter++) {
		printf("working frame %i\n", frameCounter);
		float currentFrame = float(frameCounter) / float(fps);

		constexpr float mypi = glm::pi<float>();
		vec3 eye = vec3(sin(currentFrame*2*mypi)*5.0f, 0.0f, cos(currentFrame * 2 * mypi)*5.0f);

		mat4 view = glm::lookAt(eye, vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));

		//printf("clearing buffers\n");
		clearBuffers();

		vector<Triangle> triangles;

		//printf("drawing some mesh\n");
		for (Mesh& m : backpack.children) {
			doTriangles(m, view, triangles);
			//drawMesh(m, projection, view, m.transform);
		}

		float ratio = float(frameX) / float(frameY);
		float fovY = radians(70.0f);
		float fovX = fovY * ratio;
		vec3 intoScreen(0.0f, 0.0f, -1.0f);


		/*concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;
			vec3 bary;
			vec3 boxCoord;
			vec3 hitPoint;
			//if (y == 0) { printf("on pixel %i, %i\n", x, y); }
			vec2 clipCoords = frameToNDC(ivec2(x, y));
			vec3 rayVector = glm::normalize(vec3(clipCoords * mypi, -mypi));
			Ray ray(vec3(0), rayVector);
			double hitDepth = INFINITY;
			for (const Triangle& tri : triangles) {
				if (rayHit(ray, tri.minBounding, tri.maxBounding, boxCoord)) {
					if (rayHit(ray, tri, bary, hitPoint)) {
						double distance = glm::length2(boxCoord);
						if (distance < hitDepth) {
							vec2 fragUV = bary.x * tri[0].texCoords + bary.y * tri[1].texCoords + bary.z * tri[2].texCoords;
							vec3 fragColor = tri.mesh->diffuse.sample(fragUV);
							frameBuffer[x + (y * frameX)] = fragColor;
							hitDepth = distance;
						}
					}

				}
			}

		});*/

		printf("all triangles made, about to do kdn\n");
		const KDBranch* theKDN = new KDBranch(triangles);
		printf("made kdn\n");

		concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;
			vec3 bary;
			vec3 boxCoord;
			vec3 hitPoint;
			//if (y == 0) { printf("on pixel %i, %i\n", x, y); }
			vec2 clipCoords = frameToNDC(ivec2(x, y));
			vec3 rayVector = glm::normalize(vec3(clipCoords * mypi, -mypi));
			Ray ray(vec3(0), rayVector);
			double hitDepth;
			const Triangle* tri = rayHit(theKDN, ray, bary, hitDepth);
			if (tri != NULL) {
				vec2 fragUV = bary.x * (*tri)[0].texCoords + bary.y * (*tri)[1].texCoords + bary.z * (*tri)[2].texCoords;
				vec3 fragColor = tri->mesh->diffuse.sample(fragUV);
				frameBuffer[x + (y * frameX)] = fragColor;
			}
		});
		printf("did parallelForloop\n");


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();

		saveImage((char*)("out/"+std::to_string(frameCounter) + ".png").c_str(), window);
		cin.get();

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}