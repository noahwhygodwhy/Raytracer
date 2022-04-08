#include "Raytracer.hpp"
#include "CookTorrance.hpp"
#include "clTypeDefs.hpp"

typedef float4 fvec4;
typedef float3 fvec3;

#define CPP
#include "sharedStructs.cl"

using namespace std;
using namespace std::filesystem;
using namespace glm;

#define MAX_SHAPES 100000
#define MAX_MATERIALS 10

#define MAX_PATH 100u
//#define OUTPUTPASSES 50
//#define OUTPUTFRAMES 189
//#define EVERYFRAME INFINITY
//#define CONCURRENT_FOR
//#define KDTRACE
#define CIN
//#define PPCIN

#define PIXEL_MULTISAMPLE_N 1
#define MONTE_CARLO_SAMPLES 1


//#define BASIC_BITCH


bool prd = false; //print debuging for refraction

uint32_t frameX = 1000;
uint32_t frameY = 1000;
double frameRatio = double(frameX) / double(frameY);


//fvec4 clearColor(0.21, 0.78, 0.95, 1.0);

fvec4 clearColor(0.0, 0.0, 0.0, 1.0);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame
string saveFileDirectory = "";


constexpr double bias = 1e-4;


void frameBufferSizeCallback(GLFWwindow* window, uint64_t width, uint64_t height) {
	glViewport(0, 0, GLsizei(width), GLsizei(height));
}

bool frontFacing(vec3 a, vec3 b, vec3 c) {
	mat3 m(a.x, b.x, c.x, a.y, b.y, c.y, 1.0f, 1.0f, 1.0f);
	return glm::determinant(m) > 0;
}

//random double from 0 to 1
double randDubTwo() {
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
//random double between -1 and 1
double randDubThree() {
	return (((static_cast <double> (rand()) / static_cast <double> (RAND_MAX))*2.0)-1.0);
}

dvec3 randomHemisphericalVector(dvec3 normal) {
	dvec3 outVec(0);
	normal = glm::normalize(normal);
	do {
		outVec = glm::normalize(dvec3(randDubThree(), randDubThree(), randDubThree()));
	} while (glm::dot(normal, outVec) < 0.0);
	return outVec;
}
dvec3 randomSpecularVector(dvec3 normal) {
	dvec3 outVec(0);
	normal = glm::normalize(normal);
	do {
		outVec = glm::normalize(dvec3(randDubThree(), randDubThree(), randDubThree()));
	} while (glm::dot(normal, outVec) < 0.7);
	return outVec;
}

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


AABB redoAABBs(vector<UShape>& shapes, const vector<Vertex>& vertices) {

	AABB toReturn = { fvec4(0.0), fvec4(0.0) };
	for (UShape& s : shapes) {
		switch (s.type) {
		case 0:
			fvec4 origin = fvec4(s.values.xyz, 0.0f);
			//float radius = s.values.w;
			s.boundingBox.min = origin - fvec4(fvec3(s.values.w), 0.0f);
			s.boundingBox.max = origin + fvec4(fvec3(s.values.w), 0.0f);
			break;
		case 1:
			s.boundingBox.min = s.boundingBox.max = vertices[int(s.values.x)].position;

			s.boundingBox.min = min(vertices[int(s.values.y)].position, s.boundingBox.min) -fvec4(fvec3(glm::epsilon<float>()* 10.0f), 0.0);
			s.boundingBox.min = min(vertices[int(s.values.z)].position, s.boundingBox.min) - fvec4(fvec3(glm::epsilon<float>() * 10.0f), 0.0);

			s.boundingBox.max = max(vertices[int(s.values.y)].position, s.boundingBox.max) + fvec4(fvec3(glm::epsilon<float>() * 10.0f), 0.0);
			s.boundingBox.max = max(vertices[int(s.values.z)].position, s.boundingBox.max) + fvec4(fvec3(glm::epsilon<float>() * 10.0f), 0.0);
			break;
		default:
			printf("redo aabbs bad shape type\n");
			exit(-1);
			break;
		}
		toReturn.min = glm::min(toReturn.min, s.boundingBox.min);
		toReturn.max = glm::max(toReturn.max, s.boundingBox.max);
	}

	return toReturn;
}


UShape makeSphere(fvec3 origin, float radius, uint materialIdx) {
	return UShape(fvec4(origin, radius), AABB(), materialIdx, 0u);
}
UShape makeTriangle(uint a, uint b, uint c, uint materialIdx) {
	return UShape(fvec4(a, b, c, 0), AABB(), materialIdx, 1u);
}

int main()
{

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

	cl_uint numPlatforms;
	cl_int status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetPlatformIDs failed (1)\n");
	//printf("Number of Platforms = %d\n", numPlatforms);
	cl_platform_id* platforms = new cl_platform_id[numPlatforms];
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetPlatformIDs failed (2)\n");
	cl_uint numDevices;
	cl_device_id* devices;



	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);

	devices = new cl_device_id[numDevices];
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);

	cl_device_type type;
	size_t sizes[3] = { 0, 0, 0 };
	clGetDeviceInfo(devices[0], CL_DEVICE_TYPE, sizeof(type), &type, NULL);
	
	cl_device_id device = devices[0];


	char* ui;
	size_t valueSize;
	clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(ui), &ui, NULL);

	clGetDeviceInfo(device, CL_DEVICE_VERSION, 0, NULL, &valueSize);
	ui = (char*)malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_VERSION, valueSize, ui, NULL);
	printf(" %d.%d Hardware version: %s\n", 1, 1, ui);
	free(ui);

	//printf("%c\n", ui);
	//exit(0);



	cl_context_properties properties[] = {
		CL_GL_CONTEXT_KHR, (cl_context_properties)glfwGetWGLContext(window),
		CL_WGL_HDC_KHR, (cl_context_properties)GetDC(glfwGetWin32Window(window)),
		CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[0],
	0};

	cl_context context = clCreateContext(properties, 1, &device, NULL, NULL, &status);
	cl_command_queue_properties* qProperties = new cl_command_queue_properties();
	cl_command_queue cmdQueue = clCreateCommandQueueWithProperties(context, device, qProperties, &status);

	cl_mem clOtherData = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(OtherData), NULL, &status);
	cl_mem clShapes = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(UShape) * MAX_SHAPES, NULL, &status);
	cl_mem clMaterials = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Material) * MAX_MATERIALS, NULL, &status);
	cl_mem clVerts = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Vertex) * MAX_SHAPES * 3, NULL, &status);
	cl_mem clRandomBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, frameX * frameY * sizeof(uint64_t), NULL, &status);

	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, frameX, frameY, 0, GL_RGBA, GL_FLOAT, NULL);

	cl_mem clFrameTexture = clCreateFromGLTexture(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, frameTexture, &status);
	//printf("cltesttexture status: %i\n", status);



	long unsigned int length;
	ifstream stream("Renderer.cl", ios::in | ios::ate | ios::binary);
	//stream.seekg(0, ios::end);
	length = long unsigned int(stream.tellg());
	stream.seekg(0, ios::beg);
	char* shaderSource = new char[length + 1];
	shaderSource[length] = '\0';
	const char** shaderSourcesArray = new const char* [1];
	stream.read(shaderSource, length);

	shaderSourcesArray[0] = shaderSource;

	cl_program program = clCreateProgramWithSource(context, 1, shaderSourcesArray, NULL, &status);

	const char* options = { "-cl-single-precision-constant" };
	status = clBuildProgram(program, 1, &device, options, NULL, NULL);

	if (status != CL_SUCCESS)
	{ // retrieve and print the error messages:
		size_t size;

		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size);
		cl_char* log = new cl_char[size];
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, size, log, NULL);
		printf("clBuildProgram failed:\n%s\n", log);
		exit(-1);
	}

	cl_kernel kernel = clCreateKernel(program, "render", &status);

	status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &clOtherData);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &clShapes);
	status = clSetKernelArg(kernel, 2, sizeof(cl_mem), &clVerts);
	status = clSetKernelArg(kernel, 3, sizeof(cl_mem), &clMaterials);
	status = clSetKernelArg(kernel, 4, sizeof(cl_mem), &clRandomBuffer);
	status = clSetKernelArg(kernel, 5, sizeof(cl_mem), &clFrameTexture);

	//v this is the number of items to do, so like framex and framey?
	//size_t globalWorkSize[3] = { sizes[0], sizes[1], sizes[2]};

	//values in this must be divisible by values in local work size
	//cannot be values larger than an unsigned in contained in a number of bits queried by CL_DEVICE_ADDRESS_BITS
	size_t globalWorkSize[3] = { frameX, frameY, 1 };

	//this is describing how wide the processing can be, so we set it to max the gpu can do
	//size_t localWorkSize[3] = { sizes[0], sizes[1], sizes[2]};


	//must be a 3d array volume < CL_DEVICE_MAX_WORK_GROUP_SIZE (in our case 1024)
	//and each side must be less than the coresponding size in CL_DEVICE_MAX_WORK_ITEM_SIZES (1024x1024x64)
	size_t localWorkSize[3] = { NULL, NULL, NULL };







	Shader shader("vert.glsl", "frag.glsl");
	shader.use();




	mt19937_64 numGen;
	randomBuffer = new uint64_t[frameX * frameY]();
	for (size_t i = 0; i < frameX * frameY; i++) {
		randomBuffer[i] = numGen();
		//printf("random buffer %i: %zu\n", i, randomBuffer[i]);
	}
	//exit(0);
	



	vector<Material> materials;
	materials.reserve(MAX_MATERIALS);
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 2u, 0u));//checkers
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.6f, 0.6f, 0.6f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//white light
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.54, 0.95f, 0.0f, 0.0f, 0u, 0u));//transparenty
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.9f, 0.9f, 0u, 0u));//mirrorA
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.05f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 1u));//Fog?
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(6.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//red light


	vector<UShape> shapes;
	shapes.reserve(MAX_SHAPES);

	shapes.push_back(makeSphere(fvec4(0.0f, 8.0f, 0.0f, 0.0f), 5.0f, 1u));//light

	shapes.push_back(makeSphere(fvec4(0.0f, 0.0f, 0.0f, 0.0f), 3.0f, 0u));//smal

	//spheres.push_back(Sphere(fvec4(3.0f, 3.0f, 0.0f, 0.0f), Shape(AABB(), 3u, 0), 3.0f));//transparent
	//spheres.push_back(Sphere(fvec4(0.0f, 5.0f, 5.0f, 0.0f), Shape(AABB(), 2u, 0), 3.0f));//mirror

	//spheres.push_back(Sphere(fvec4(0.0f, 3.0f, 0.0f, 0.0f), Shape(AABB(), 4u, 0), 3.0f));//foggy sphere
	//spheres.push_back(Sphere(fvec4(3.6f, 1.5f, 0.0f, 0.0f), Shape(AABB(), 5u, 0), 0.5f));//red light




	//vector<Triangle> triangles;
	//triangles.reserve(MAX_TRIANGLES);

	vector<Vertex> vertices;
	vertices.reserve(MAX_SHAPES * 3);

	/*
	vertices.push_back(Vertex(fvec4(-10, 0, -15, 0), fvec4(0, 1, 0, 0), vec4(0, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(10, 0, -15, 0), fvec4(0, 1, 0, 0), vec4(1, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(10, 0, 15, 0), fvec4(0, 1, 0, 0), vec4(1, 1, 0, 0)));
	vertices.push_back(Vertex(fvec4(-10, 0, 15, 0), fvec4(0, 1, 0, 0), vec4(0, 1, 0, 0)));*/

	vertices.push_back(Vertex(fvec4(-5, 0, -15, 0), fvec4(0, 1, 0, 0), vec4(0, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(15, 0, -15, 0), fvec4(0, 1, 0, 0), vec4(1, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(15, 0, 15, 0), fvec4(0, 1, 0, 0), vec4(1, 1, 0, 0)));
	vertices.push_back(Vertex(fvec4(-5, 0, 15, 0), fvec4(0, 1, 0, 0), vec4(0, 1, 0, 0)));

	shapes.push_back(makeTriangle(0, 3, 2, 0u));
	shapes.push_back(makeTriangle(0, 2, 1, 0u));

	//triangles.push_back(Triangle(Shape(AABB(), 0u, 1), 0, 3, 2));
	//triangles.push_back(Triangle(Shape(AABB(), 0u, 1), 0, 2, 1));






	uint32_t frameCounter = 0;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	uint32_t fps = 30;


	cl_event* waitAfterWrites = new cl_event[3];

	AABB sceneBounding = redoAABBs(shapes, vertices);




	status = clEnqueueWriteBuffer(cmdQueue, clShapes, CL_FALSE, 0, sizeof(Shape) * MAX_SHAPES, shapes.data(), 0, NULL, &waitAfterWrites[0]);
	status = clEnqueueWriteBuffer(cmdQueue, clVerts, CL_FALSE, 0, sizeof(Vertex) * MAX_SHAPES * 3, vertices.data(), 0, NULL, &waitAfterWrites[1]);
	status = clEnqueueWriteBuffer(cmdQueue, clRandomBuffer, CL_FALSE, 0, sizeof(uint64_t) * frameX * frameY, randomBuffer, 0, NULL, &waitAfterWrites[2]);
	cl_event* waitAfterFinalWrite = new cl_event;
	status = clEnqueueWriteBuffer(cmdQueue, clMaterials, CL_TRUE, 0, sizeof(Material) * MAX_MATERIALS, materials.data(), 3, waitAfterWrites, waitAfterFinalWrite);

	cl_event* otherDataEvent = new cl_event;
	cl_event* waitAfterProcessing = new cl_event;

#if defined(OUTPUTFRAMES)|| defined(OUTPUTPASSES)

	time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/"+saveFileDirectory).c_str());
	cl_event* readEvent;


#ifdef OUTPUTFRAMES
	for (; frameCounter < OUTPUTFRAMES;) {
	double currentFrame = double(frameCounter) / double(fps);



#endif
#ifdef OUTPUTPASSES
		{
#endif
#else

	while (!glfwWindowShouldClose(window)) {

#ifdef EVERYFRAME

		double currentFrame = double(frameCounter) / double(fps);
#else
		double currentFrame = glfwGetTime();
#endif

#endif
		frameCounter++;



		//currentTime = currentFrame;
		deltaTime = currentFrame - lastFrame;
		printf("that frame took %f seconds\n", deltaTime);
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


		//fvec3 eye = fvec3(sin(currentFrame) * 16, 4, cos(currentFrame) * 16);

		fvec3 eye = fvec3(0.0f, 5.0f, 6.0f);
		fvec3 lookat = fvec3(0.0, 0.0, 0.0);

		fvec3 camForward = glm::normalize(lookat - eye);
		fvec3 camUp = glm::normalize(fvec3(0.0, 1, 0.0));
		fvec3 camRight = glm::cross(camForward, camUp);
		camUp = glm::cross(camRight, camForward);

		float viewPortHeight = 2.0f;
		float viewPortWidth = viewPortHeight * frameRatio;

		float fov = 120;
		float focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));


		sceneBounding = redoAABBs(shapes, vertices);
		OtherData otherData = { 
			clearColor, 
			fvec4(eye, 0.0f), 
			fvec4(camRight, 0.0f),
			fvec4(camUp, 0.0f),
			fvec4(camForward, 0.0),
			focal,
			currentFrame,
			MAX_PATH,
			uint(shapes.size()),
			MONTE_CARLO_SAMPLES
		};

		status = clEnqueueWriteBuffer(cmdQueue, clOtherData, CL_FALSE, 0, sizeof(OtherData), &otherData, 1, waitAfterFinalWrite, otherDataEvent);
		status = clEnqueueNDRangeKernel(cmdQueue, kernel, 2, NULL, globalWorkSize, NULL, 1, otherDataEvent, waitAfterProcessing);

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


