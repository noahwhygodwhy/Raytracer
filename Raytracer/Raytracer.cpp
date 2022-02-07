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
#include "DirectionalLight.hpp"


using namespace std;
using namespace glm;

#define MAXBOUNCES 8
//#define OUTPUTFRAMES

uint32_t framesToRender = 151;


uint32_t frameX = 500;
uint32_t frameY = 500;
double frameRatio = double(frameX) / double(frameY);

//dvec2 pixelOffset = vec2(0.5 / float(frameX), 0.5 / float(frameY));

float screenRatio = float(frameX) / float(frameY);
uint64_t screenY = glm::max(1000u, frameX);// scaleFactor* frameY;
uint64_t screenX = uint64_t(glm::max(1000u, frameX) * screenRatio);// scaleFactor* frameX;

//const float near = 0.1f;
//const float far = 50.0f;

dvec3 clearColor(0.21, 0.78, 0.95);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame


struct FrameInfo {
	vector<Shape*> shapes;
	vector<Light*> lights;
	dmat4 view;
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
					minRayResult = rayResult;
				}
			}
		//}
	}
	return minRayResult;
}




bool rayTrace(const Ray& ray, const FrameInfo& fi, dvec3& colorResult, uint32_t layer = 0)
{
	double bias = 1e-4;
	colorResult = dvec3(0);
	HitResult minRayResult = shootRay(ray, fi);
	if (minRayResult.shape != NULL) {
		dvec3 diffuse = dvec3(0.0);
		dvec3 specular = dvec3(0.0);
		dvec3 ambient = (clearColor*0.1)*minRayResult.shape->mat.ka;
		for (Light* light : fi.lights) {

			Ray lightRay = light->getRay(minRayResult.position+(minRayResult.normal*bias), fi.view);
			HitResult lightRayResult = shootRay(lightRay, fi);

			double lightDistance = light->getDistance(minRayResult.position, fi.view);
			if (lightDistance < lightRayResult.depth) {//shawows

				double attenuation = light->getAttenuation(lightDistance);

				dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal)- lightRay.direction);

				//ambient += minRayResult.shape->mat.ka;

				double spec = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, ray.inverseDirection)), minRayResult.shape->mat.ns);/*^to the specular exponent*/
				specular += (minRayResult.shape->mat.ks*light->color * spec) / attenuation;

				double diff = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
				diffuse += (minRayResult.shape->mat.kd*light->color * diff) / attenuation;

			}
		}

		bool oddSecond = uint32_t(floor(fi.currentTime)) % 2 == 0;

		colorResult = (ambient + diffuse + specular) * minRayResult.shape->getColor(minRayResult);
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

	//shapes.push_back(new Sphere(dvec3(0), 1.0, oscilateX));
	



	//TODO: two spheres makes it hard. NEED to do the octtree

	shapes.push_back(new Sphere(dvec3(0, 3, 0), 2.0, materials.at("Bronze"), oscilateX));
	shapes.push_back(new Sphere(dvec3(0, -3, 0), 2.0, materials.at("Polished Bronze"), oscilateX));
	//shapes.push_back(new Sphere(dvec3(3, 0, 0), 1.0, noMovement));

	vector<Light*> lights;

	lights.push_back(new PointLight(
		dvec3(0.0, 7.0, 7.0),
		dvec3(0.8, 0.8, 1.0),
		dvec3(1.0, 0.09, 0.032)
	));

	/*lights.push_back(new DirectionalLight(
		dvec3(-0.3, -1.0, 0.0),
		dvec3(0.5, 0.0, 0.0)
	));*/




	FrameInfo fi;
	fi.currentTime = 0;
	fi.shapes = shapes;
	fi.lights = lights;
	fi.view = mat4(1.0);



	uint32_t fps = 24;
#ifdef OUTPUTFRAMES
	for (uint32_t frameCounter = 0; frameCounter < framesToRender; frameCounter++) {
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
		//vec3 eye = vec3(sin(currentFrame) * 10.0f, 2.0f, cos(currentFrame) * 10.0f);

		clearBuffers();


		dvec3 cameraForward = dvec3(0.0, 0.0, -1.0);
		dvec3 cameraUp = dvec3(0.0, 1.0, 0.0);

		dvec3 eye = vec3(0.0, 2.0, 10.0); 
		//dmat4 view = glm::lookAt(eye, eye+cameraForward, cameraUp);
		fi.view = glm::lookAt(eye, eye + cameraForward, cameraUp);
		//printf("focal: %f\n", focal);



		double pixelWidth = 2.0 / double(frameX);
		double halfPixel = 1.0 / double(frameX);
		double quarterPixel = 0.5 / double(frameX);

		double startingX = (-pixelWidth * (frameX / 2.0f));

		double viewPortHeight = 2.0f;
		double viewPortWidth = viewPortHeight * frameRatio;

		double fov = 70;
		//printf("current frame: %f\n", fov);
		double focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));
		//double focal = 1.0;


		//printf("Framerati: %f\n", frameRatio);

		dvec3 horizontal = dvec3(viewPortWidth, 0, 0);
		dvec3 vertical = dvec3(0, viewPortHeight, 0);

		dvec2 hv(viewPortWidth, viewPortHeight);
		dvec3 origin(0);
		dvec3 lowerLeftCorner = origin - (horizontal / 2.0) - (vertical / 2.0) - dvec3(0, 0, focal);



		//printf("lowerleftcorner: %s\n", glm::to_string(lowerLeftCorner).c_str());



		//for (uint64_t i = 0; i < frameX * frameY; i++) {
		concurrency::parallel_for(uint64_t(0), uint64_t(frameX * frameY), [&](uint64_t i) {
			uint32_t x = i % frameX;
			uint32_t y = i / frameX;

			//double u = double(x) / double((frameX - 1.0));
			//double v = double(y) / double((frameY - 1.0));
			
			dvec2 uv(double(x) / double((frameX - 1.0)), double(y) / double((frameY - 1.0)));
			dvec2 clipSpacePixelSize(dvec2(1.0 / double(frameX - 1.0), 1.0 / double(frameY - 1.0)));



			//dvec3 rayVector = lowerLeftCorner + (uv.x * horizontal) + (uv.y * vertical);

			//Ray initialRay(dvec3(0.0, 0.0, focal), glm::normalize(rayVector));
			//dvec3 colorResult(0.0);


			//multisample n*n
			uint32_t n = 1;
			dvec3 colorAcum(0);

			dvec2 pixelBottomLeft = ((uv.x * horizontal) + (uv.y * vertical)).xy - (clipSpacePixelSize / 2.0);

			for (uint32_t x = 1; x <= n; x++) {
				double offsetX = x * (clipSpacePixelSize.x / (n + 1));
				for (uint32_t y = 1; y <= n; y++) {
					double offsetY = y * (clipSpacePixelSize.y / (n + 1));
					//printf("x: %u, y: %u\n", x, y);
					//printf("offsetY: %f, offsetX: %f\n", offsetX, offsetY);
					dvec3 rayVector = lowerLeftCorner + (uv.x * horizontal) + (uv.y * vertical) + dvec3(offsetX, offsetY, 0.0);

					dvec3 colorOut;

					Ray initialRay(dvec3(0.0, 0.0, focal), glm::normalize(rayVector));
					bool hit = rayTrace(initialRay, fi, colorOut);
					if (hit) {
						colorAcum += colorOut;
					}
					else {

						colorAcum += clearColor;
					}

				}
			}
			frameBuffer[x + (y * frameX)] = colorAcum / double(n * n);
		});
		//}

		//printf("did parallelForloop\n");


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();
#ifdef OUTPUTFRAMES
		saveImage((char*)("out/"+std::to_string(frameCounter) + ".png").c_str(), window);
#endif
		//cin.get();

	}
	//pool.stop();
	std::printf("closing\n");
	glfwTerminate();
}


