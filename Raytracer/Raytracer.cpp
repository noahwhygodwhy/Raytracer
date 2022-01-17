#include <iostream>
#include <vector>
#include <fstream>

//#define GLM_FORCE_LEFT_HANDED 
//#define GLM_RIGHT_HANDED 

#include "glad/glad.h"
#include <GLFW/glfw3.h>
#define GLM_SWIZZLE 
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "glm/gtx/string_cast.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION  
#include "stb_image_write.h"

#include "Shader.hpp"

using namespace std;
using namespace glm;



const uint64_t screenX = 1000;
const uint64_t screenY = 1000;
const fvec2 pixelOffset = vec2(0.5 / float(screenX), 0.5 / float(screenY));


const float near = 0.1f;
const float far = 10.0f;

vec3 clearColor(0.0f, 0.0f, 0.0f);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame

/*struct Vertex {
	vec3 position;
	vec3 normal;
	vec3 color;
};*/

struct Point {
	fvec3 position;
	fvec3 color;
};

struct Triangle {
	Vertex vertices[3];
	ivec2 maxFrags;
	ivec2 minFrags;
	ivec2 zLimits;
	bool onscreen = true;

	vec2 v0;//for barycentry coords
	vec2 v1;//^
	float d00;
	float d01;
	float d11;
	float denom;
	Triangle(Vertex a, Vertex b, Vertex c) {
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;

		v0 = b.position.xy - a.position.xy;
		v1 = c.position.xy - a.position.xy;

		d00 = glm::dot(v0, v0);
		d01 = glm::dot(v0, v1);
		d11 = glm::dot(v1, v1);

		denom = (d00 * d11) - (d01 * d01);


		/*printf("a: %s\n", glm::to_string(a.position).c_str());
		printf("b: %s\n", glm::to_string(b.position).c_str());
		printf("c: %s\n", glm::to_string(c.position).c_str());*/

		//maxFrags.x = uint32_t(ceil(glm::max(glm::max(a.position.x, b.position.x), c.position.x)));

		maxFrags.x = uint32_t(ceil(glm::min((glm::max(glm::max(a.position.x, b.position.x), c.position.x)+1.0f)/2.0f * screenX, float(screenX-1))));
		minFrags.x = uint32_t(floor(glm::max((glm::min(glm::min(a.position.x, b.position.x), c.position.x) + 1.0f) / 2.0f * screenX, 0.0f)));
		maxFrags.y = uint32_t(ceil(glm::min((glm::max(glm::max(a.position.y, b.position.y), c.position.y) + 1.0f) / 2.0f * screenY, float(screenY - 1))));
		minFrags.y = uint32_t(floor(glm::max((glm::min(glm::min(a.position.y, b.position.y), c.position.y) + 1.0f) / 2.0f * screenY, 0.0f)));

		zLimits.x = glm::min(glm::min(a.position.z, b.position.z), c.position.z);
		zLimits.y = glm::max(glm::max(a.position.z, b.position.z), c.position.z);

		if (minFrags.x > screenX || maxFrags.x < 0 || minFrags.y > screenY || maxFrags.y < 0) { // || zLimits.x > far || zLimits.y < near) { //TODO: get this working, idk why it's wrong
			onscreen = false;
		}
	}
	Vertex operator[](uint32_t index) const {
		return vertices[index];
	}

	vec3 getBaryCoords(vec2 p) const {
		vec2 v2 = p - vertices[0].position.xy;
		float d20 = dot(v2, v0);
		float d21 = dot(v2, v1);

		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		return vec3(1.0f - v - w, v, w);
	}


};


float randFloat() {
	return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
}


vector<Vertex> getIcosahedron()
{

	float r = 0.5f;
	float goldenRatio = (1.0f + glm::sqrt(5.0f)) / 2.0f;
	float sl = r / goldenRatio;

	//points on x plane
	vec3  x1 = vec3(0, r, sl);
	vec3 x2 = vec3(0, -r, sl);
	vec3 x3 = vec3(0, -r, -sl);
	vec3 x4 = vec3(0, r, -sl);
	//x4 = [0,0,0]

	//points on y plane

	vec3 y1 = vec3(sl, 0, r);
	vec3 y2 = vec3(sl, 0, -r);
	vec3 y3 = vec3(-sl, 0, -r);
	vec3 y4 = vec3(-sl, 0, r);

	//points on z plane

	vec3 z1 = vec3(r, sl, 0);
	vec3 z2 = vec3(-r, sl, 0);
	vec3 z3 = vec3(-r, -sl, 0);
	vec3 z4 = vec3(r, -sl, 0);

	vector<vec3> orderedPoints = {
		x1, x4, z2,
		x4, x1, z1,
		x2, x3, z4,
		x3, x2, z3,
		y1, y4, x2,
		y4, y1, x1,
		y2, y3, x4,
		y3, y2, x3,
		z1, z4, y2,
		z4, z1, y1,
		z2, z3, y4,
		z3, z2, y3,
		x1, y1, z1,
		x4, z1, y2,
		x4, y3, z2,
		x1, z2, y4,
		x2, z4, y1,
		x2, y4, z3,
		x3, y2, z4,
		x3, z3, y3
	};

	//the normal is the point lol
	//the color is based on it's sin(x), cos(z), tan(y)* just for fun
	vector<Vertex> toReturn;
	vec3 red(1.0f, 0.0f, 0.0f);
	//for (const vec3& p : orderedPoints) {//does sphere normals
	//	toReturn.push_back({ p, p, red/*vec3(p.x + 0.5f, p.z + 0.5f, p.y + 0.5f)*/ });
	//}
	for (uint32_t i = 0; i < orderedPoints.size(); i += 3) { //does sharp ico normals;

		vec3 color(randFloat(), randFloat(), randFloat());
		vec3 centroid = (orderedPoints[i+0] + orderedPoints[i+1] + orderedPoints[i + 2]) / 3.0f;
		vec3 centroidNormal = glm::normalize(centroid);
		toReturn.push_back({ orderedPoints[i + 0], centroidNormal, red/*vec3(p.x + 0.5f, p.z + 0.5f, p.y + 0.5f)*/});
		toReturn.push_back({ orderedPoints[i + 1], centroidNormal, red/*vec3(p.x + 0.5f, p.z + 0.5f, p.y + 0.5f)*/ });
		toReturn.push_back({ orderedPoints[i + 2], centroidNormal, red/*vec3(p.x + 0.5f, p.z + 0.5f, p.y + 0.5f)*/ });
	}
	return toReturn;
}



//TODO: do i need to transform the normals as well..?
//Applies the projection and view matrix to the point
Vertex transformIt(const mat4& projection, const mat4& view, const mat4& model, Vertex v) {


	mat4 mvp = projection * view * model;
	vec4 newPos = mvp* vec4(v.position, 1.0f);
	v.position = vec3(newPos.x, newPos.y, newPos.z)/newPos.w;
	return v;

	if (uint32_t(lastFrame) % 2 == 0) {

		v.position = projection * vec4(v.position, 1.0f);
	}
	else {
		v.position = projection * view * vec4(v.position, 1.0f);
	}

}




//sign() and pointInTriangle() based off of the code in this stackoverflow post https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float sign(vec2 a, vec2 b, vec2 c) {
	return ((a.x - c.x) * (b.y - c.y)) - ((b.x - c.x) * (a.y - c.y));
}

bool pointInTriangle(vec2 point, const Triangle& tri) {

	float signab, signbc, signca;
	bool hasNeg, hasPos;

	signab = sign(point, tri[0].position, tri[1].position);
	signbc = sign(point, tri[1].position, tri[2].position);
	signca = sign(point, tri[2].position, tri[0].position);
	
	hasNeg = (signab < 0) || (signbc < 0) || (signca < 0);
	hasPos = (signab > 0) || (signbc > 0) || (signca > 0);

	return !(hasNeg && hasPos);
}





void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

}



float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}

int main()
{

	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SRGB_CAPABLE, 1);
	glfwWindowHint(GLFW_SAMPLES, 16);
	GLFWwindow* window = glfwCreateWindow(screenX, screenY, "Renderer", NULL, NULL);
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
	glViewport(0, 0, screenX, screenY);
	glEnable(GL_FRAMEBUFFER_SRGB);


	unsigned int VBO, VAO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);

	glClearColor(clearColor.x, clearColor.y, clearColor.z, 1.0f);





	Shader shader("vert.glsl", "frag.glsl");
	shader.use();

	vector<Vertex> vertices = getIcosahedron();/* {
		{vec3(-1.0f, -0.1f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.0f, 0.0f)},
		{vec3(0.2f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f)},
		{vec3(1.0f, 0.0f, 0.2f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f)},
	};*/

	//make framebuffers
	vec3* frameBuffer = new vec3[screenX * screenY]();
	float* depthBuffer = new float[screenX * screenY]();

	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	while (!glfwWindowShouldClose(window))
	{

		//frametime calculation
		double currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		//projection/view matrix 
		float ratio = float(screenX) / float(screenY);
		mat4 projection = glm::perspective(radians(70.0f), ratio, near, far);
		vec3 eye = vec3(sin(currentFrame)*2.0f, 0.0f, cos(currentFrame)*2.0f);
		mat4 view = glm::lookAt(eye, vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));
		mat4 model = glm::rotate(mat4(1.0f), sinf(currentFrame), vec3(1.0f, 0.0f, 0.0f));

		//transform vertices to account for projection, view, and model
		vector<Vertex> transformedVertices;
		transformedVertices.resize(vertices.size());
		for (uint64_t i = 0; i < vertices.size(); i++) {
			transformedVertices[i] = transformIt(projection, view, model, vertices[i]);
		}

		//construct triangles out of vertices
		vector<Triangle> triangles;
		for (int i = 0; i < transformedVertices.size(); i += 3) {
			triangles.push_back(Triangle(
				transformedVertices[i + 0],
				transformedVertices[i + 1],
				transformedVertices[i + 2]
			));
		}
		transformedVertices.clear();
		transformedVertices.shrink_to_fit();

		//clear the color and depth buffers from last frame
		for (uint64_t x = 0; x < screenX; x++) {
			for (uint64_t y = 0; y < screenY; y++) {
				frameBuffer[x + (y * screenX)] = clearColor;
				depthBuffer[x + (y * screenX)] = INFINITY;
			}
		}


		vec3 lightDirection(-0.3f, -1.0f, -0.2f);
		vec3 lightColor(0.8f, 0.8f, 0.8f);
		vec3 ambientLight(0.2f, 0.2f, 0.2f);
		vec3 specular;
		vec3 reversedLightDirection = glm::normalize(-lightDirection);

		for (const Triangle& tri : triangles) {
			if (tri.onscreen) {
				for (uint32_t fragx = tri.minFrags.x; fragx <= tri.maxFrags.x; fragx++) { //Lines that try to limit the box of fragments to where the triangle should be
					for (uint32_t fragy = tri.minFrags.y; fragy <= tri.maxFrags.y; fragy++) {
				//for (uint32_t fragx = 0; fragx < screenX; fragx++) {
					//for (uint32_t fragy = 0; fragy < screenY; fragy++) {

						fvec2 fragCoord = fvec2((floatDiv(fragx, screenX) * 2.0f) - 1.0f, (floatDiv(fragy, screenY) * 2.0f) - 1.0f) + pixelOffset;
						if (pointInTriangle(fragCoord, tri)) {

							vec3 beryCoords = tri.getBaryCoords(fragCoord);

							float fragDepth = glm::dot(vec3(tri[0].position.z, tri[1].position.z, tri[2].position.z), beryCoords);
							if (fragDepth < depthBuffer[fragx + (fragy * screenX)]) {
								vec3 fragColor = mat3(tri[0].color, tri[1].color, tri[2].color) * beryCoords;
								vec3 fragNormal = mat3(tri[0].normal, tri[1].normal, tri[2].normal) * beryCoords;

								vec3 ambientResult = fragColor * ambientLight;

								float diff = glm::max(glm::dot(fragNormal, reversedLightDirection), 0.0f);
								vec3 diffuse = diff * lightColor;
								vec3 diffuseResult = diffuse * fragColor;

								depthBuffer[fragx + (fragy * screenX)] = fragDepth;
								frameBuffer[fragx + (fragy * screenX)] = diffuseResult + ambientResult;
							}
						}
					}
				}
			}
		}

		/*unsigned char* normalizedFrameData = new unsigned char[screenX * screenY * 3];
		for (uint64_t i = 0; i < screenX * screenY; i++) {
			uint64_t nfdIndex = i * 3;
			normalizedFrameData[nfdIndex + 0] = unsigned char(floor(frameBuffer[i].x + 0.5f)) * 255;
			normalizedFrameData[nfdIndex + 1] = unsigned char(floor(frameBuffer[i].y + 0.5f)) * 255;
			normalizedFrameData[nfdIndex + 2] = unsigned char(floor(frameBuffer[i].z + 0.5f)) * 255;
		}*/



		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, screenX, screenY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		processInput(window);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	std::printf("closing\n");
	glfwTerminate();
}