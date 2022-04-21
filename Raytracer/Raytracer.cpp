#include "Raytracer.hpp"
#include "clTypeDefs.hpp"

typedef float4 fvec4;
typedef float3 fvec3;

#define CPP
#include "sharedStructs.cl"

using namespace std;
using namespace std::filesystem;
using namespace glm;

#define MAX_SHAPES 1000
#define MAX_MATERIALS 10

#define MAX_PATH 200u
#define OUTPUTFRAMES 1024
#define EVERYFRAME INFINITY
//#define CONCURRENT_FOR
//#define KDTRACE
//#define CIN


#define PIXEL_MULTISAMPLE_N 1
#define MONTE_CARLO_SAMPLES 200


//#define BASIC_BITCH


bool prd = false; //print debuging for refraction

uint32_t frameX = 400;
uint32_t frameY = 400;
double frameRatio = double(frameX) / double(frameY);


fvec4 clearColor(0.21, 0.78, 0.95, 1.0);

//fvec4 clearColor(0.0, 0.0, 0.0, 1.0);

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
	GLsizei nrChannels = 4;
	GLsizei stride = nrChannels * width;
	stride += (stride % 4) ? (4 - stride % 4) : 0;
	GLsizei bufferSize = stride * height;
	std::vector<float> buffer(bufferSize);
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer.data());


	//TODO read it as floats, apply correction, then convert it to GL_UNSIGNED_BYTE and output image
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



//INITIALIZE OPENGL STUFF
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
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


//initialize opencl devices
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

	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
	devices = new cl_device_id[numDevices];
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
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


//Set up opencl memory
	cl_context_properties properties[] = {
		CL_GL_CONTEXT_KHR, (cl_context_properties)glfwGetWGLContext(window),
		CL_WGL_HDC_KHR, (cl_context_properties)GetDC(glfwGetWin32Window(window)),
		CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[0],
	0};

	cl_context context = clCreateContext(properties, 1, &device, NULL, NULL, &status);
	printf("context status: %i\n", status);
	cl_command_queue_properties* qProperties = new cl_command_queue_properties();
	cl_command_queue cmdQueue = clCreateCommandQueueWithProperties(context, device, qProperties, &status);
	printf("cmdqueue status: %i\n", status);

	cl_mem clOtherData = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(OtherData), NULL, &status);
	printf("create buffer 0 status: %i\n", status);
	cl_mem clShapes = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(UShape) * MAX_SHAPES, NULL, &status);
	printf("create buffer 1 status: %i\n", status);
	cl_mem clMaterials = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Material) * MAX_MATERIALS, NULL, &status);
	printf("create buffer 2 status: %i\n", status);
	cl_mem clVerts = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Vertex) * MAX_SHAPES * 3, NULL, &status);
	printf("create buffer 3 status: %i\n", status);
	cl_mem clRandomBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, frameX * frameY * sizeof(uint64_t), NULL, &status);
	printf("create buffer 4 status: %i\n", status);
	cl_mem clToneMapData = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ToneMapStruct), NULL, &status);
	printf("create buffer 4 status: %i\n", status);





//Set up opengl buffer to be drawn to by opencl
	GLuint frameFBO;
	glGenFramebuffers(1, &frameFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, frameFBO);

	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, frameX, frameY, 0, GL_RGBA, GL_FLOAT, NULL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, frameTexture, 0);

	cl_mem clFrameTexture = clCreateFromGLTexture(context, CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, frameTexture, &status);
	printf("create buffer 5 status: %i\n", status);
	clEnqueueAcquireGLObjects(cmdQueue, 1, &clFrameTexture, 0, NULL, NULL);
	//printf("cltesttexture status: %i\n", status);



//read in opencl program and compile
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




//create the kernels and set memory arguments
	cl_kernel raytraceKernel = clCreateKernel(program, "render", &status);
	printf("kernal status: %i\n", status);

	status = clSetKernelArg(raytraceKernel, 0, sizeof(cl_mem), &clOtherData);
	printf("set arg 0 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 1, sizeof(cl_mem), &clShapes);
	printf("set arg 1 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 2, sizeof(cl_mem), &clVerts);
	printf("set arg 2 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 3, sizeof(cl_mem), &clMaterials);
	printf("set arg 3 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 4, sizeof(cl_mem), &clRandomBuffer);
	printf("set arg 4 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 5, sizeof(cl_mem), &clFrameTexture);
	printf("set arg 5 status: %i\n", status);
	status = clSetKernelArg(raytraceKernel, 6, sizeof(cl_mem), &clToneMapData);
	printf("set arg 6 status: %i\n", status);
	 


	cl_kernel tonemapKernel = clCreateKernel(program, "toneMap", &status);
	printf("kernal status: %i\n", status);

	status = clSetKernelArg(tonemapKernel, 0, sizeof(cl_mem), &clOtherData);
	printf("set arg2 0 status: %i\n", status);
	status = clSetKernelArg(tonemapKernel, 1, sizeof(cl_mem), &clFrameTexture);
	printf("set arg2 1 status: %i\n ", status);
	status = clSetKernelArg(tonemapKernel, 2, sizeof(cl_mem), &clFrameTexture);
	printf("set arg2 2 status: %i\n ", status);
	status = clSetKernelArg(tonemapKernel, 3, sizeof(cl_mem), &clToneMapData);
	printf("set arg2 3 status: %i\n", status);





//Documentation todo:
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




//opengl shader
	Shader shader("vert.glsl", "frag.glsl");
	shader.use();


//set up initial random buffer state cause gpu random is hard
	mt19937_64 numGen;
	randomBuffer = new uint64_t[frameX * frameY]();
	for (size_t i = 0; i < frameX * frameY; i++) {
		randomBuffer[i] = numGen();
	}
	



	vector<Material> materials;
	materials.reserve(MAX_MATERIALS);
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 2u, 0u));//checkers
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(5.0f, 5.0f, 5.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//white light
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.54, 0.95f, 0.0f, 0.0f, 0u, 0u));//transparenty
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.9f, 0.9f, 0u, 0u));//mirrorA
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0030f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 3u));//Fog?
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(6.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//red 
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 6.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//blue light
	materials.push_back(Material(fvec4(clearColor.xyz, 0.0f), fvec4(fvec3(clearColor.xyz)/*5.0f*/, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0u, 0u));//environ light


	vector<UShape> shapes;
	shapes.reserve(MAX_SHAPES);

	//FOG SHAPES
	//shapes.push_back(makeSphere(fvec4(6.7f, 0.0f, 0.0f, 0.0f), 1.0f, 5u));//light
	//shapes.push_back(makeSphere(fvec4(-5.0f, 0.0f, 0.0f, 0.0f), 1.0f, 6u));//light
	//shapes.push_back(makeSphere(fvec4(0.0f, 0.0f, 0.0f, 0.0f), 19.9f, 4u));//fog
	//shapes.push_back(makeSphere(fvec4(-12.0f, 0.0f, 0.0f, 0.0f), 3.0f, 0u));//small ball




	//shapes.push_back(makeSphere(fvec4(0.0, 7.0f, 0.0f, 0.0f), 5.0f, 4u));//fog ball
	shapes.push_back(makeSphere(fvec4(0.0f, 6.0f, 0.0f, 0.0f), 5.0f, 4u));//fog ball



	//shapes.push_back(makeSphere(fvec4(0.0, 5.0f, 0.0f, 0.0f), 5.0f, 3u));//reflective ball
	//shapes.push_back(makeSphere(fvec4(-7.5f, 6.0f, 7.5f, 0.0f), 5.0f, 2u));//refractive ball
	//shapes.push_back(makeSphere(fvec4(7.5f, 6.0f, 7.5f, 0.0f), 5.0f, 4u));//fog ball

	//shapes.push_back(makeSphere(fvec4(0.0f, 6.0f, 0.0f, 0.0f), 5.0f, 4u));//fog ball

	//shapes.push_back(makeSphere(fvec4(0.0, 7.0f, 4.0f, 0.0f), 5.0f, 2u));//reflective ball
	//shapes.push_back(makeSphere(fvec4(8.0, 7.0f, 4.0f, 0.0f), 5.0f, 4u));//fog ball
	//shapes.push_back(makeSphere(fvec4(-8.0, 5.0f, 0.0f, 0.0f), 5.0f, 3u));//refractive ball

	shapes.push_back(makeSphere(fvec4(0.0, 25.0f, 0.0f, 0.0f), 10.0f, 1u));//light ball
	shapes.push_back(makeSphere(fvec4(0.0, 0.0f, 0.0f, 0.0f), 200.0f, 7u));//environment ball





	vector<Vertex> vertices;
	vertices.reserve(MAX_SHAPES * 3);


	vertices.push_back(Vertex(fvec4(-15, 0, -30, 0), fvec4(0, 1, 0, 0), vec4(0, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(30, 0, -30, 0), fvec4(0, 1, 0, 0), vec4(1, 0, 0, 0)));
	vertices.push_back(Vertex(fvec4(30, 0, 30, 0), fvec4(0, 1, 0, 0), vec4(1, 1, 0, 0)));
	vertices.push_back(Vertex(fvec4(-15, 0, 30, 0), fvec4(0, 1, 0, 0), vec4(0, 1, 0, 0)));

	shapes.push_back(makeTriangle(0, 3, 2, 0u));
	shapes.push_back(makeTriangle(0, 2, 1, 0u));



	uint32_t frameCounter = 0;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	uint32_t fps = 12;

	cl_event* waitAfterWrites = new cl_event[5];

	AABB sceneBounding = redoAABBs(shapes, vertices);



	OtherData otherData = {
		clearColor,
		fps,
		MAX_PATH,
		uint(shapes.size()),
		MONTE_CARLO_SAMPLES
	};

	ToneMapStruct toneMapData = {
		{10000, 10000, 10000}, {0, 0, 0}, {0, 0, 0}

	};

	//write all the geometry data to gpu buffers
	status = clEnqueueWriteBuffer(cmdQueue, clShapes, CL_FALSE, 0, sizeof(UShape) * MAX_SHAPES, shapes.data(), 0, NULL, &waitAfterWrites[0]);
	printf("write 0 status: %i\n", status);
	status = clEnqueueWriteBuffer(cmdQueue, clVerts, CL_FALSE, 0, sizeof(Vertex) * MAX_SHAPES * 3, vertices.data(), 0, NULL, &waitAfterWrites[1]);
	printf("write 1 status: %i\n", status);
	status = clEnqueueWriteBuffer(cmdQueue, clRandomBuffer, CL_FALSE, 0, sizeof(uint64_t) * frameX * frameY, randomBuffer, 0, NULL, &waitAfterWrites[2]);
	printf("write 2 status: %i\n", status);
	status = clEnqueueWriteBuffer(cmdQueue, clMaterials, CL_FALSE, 0, sizeof(Material) * MAX_MATERIALS, materials.data(), 0, NULL, &waitAfterWrites[3]);
	printf("write 3 status: %i\n", status);
	status = clEnqueueWriteBuffer(cmdQueue, clOtherData, CL_FALSE, 0, sizeof(OtherData), &otherData, 0, NULL, &waitAfterWrites[4]);
	printf("write 4 status: %i\n", status);





	clWaitForEvents(5, waitAfterWrites);




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
		//frame time and fps calculation
		frameCounter++;
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





		//set the frame counter argument
		status = clSetKernelArg(raytraceKernel, 7, sizeof(cl_uint), &frameCounter);
		printf("set arg 7 status: %i\n", status);
		
		status = clEnqueueWriteBuffer(cmdQueue, clToneMapData, CL_TRUE, 0, sizeof(ToneMapStruct), &toneMapData, 0, NULL,NULL);
		printf("write tonemap status: %i\n", status);


		//do the raytracing pass on gpu
		cl_event* firstPassEvent = new cl_event();
		status = clEnqueueNDRangeKernel(cmdQueue, raytraceKernel, 2, NULL, globalWorkSize, NULL, 0, NULL, firstPassEvent);
		printf("range kernel: %i\n", status);
			
		//do the tonemapping pass on gpu
		/*cl_event* toneMapEvent = new cl_event();
		status = clEnqueueNDRangeKernel(cmdQueue, tonemapKernel, 2, NULL, globalWorkSize, NULL, 1, firstPassEvent, toneMapEvent);
		printf("range kernel: %i\n", status);

		clWaitForEvents(1, toneMapEvent);*/


		//draw it to the glfw window
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glDrawBuffer(GL_BACK);

		glBindFramebuffer(GL_READ_FRAMEBUFFER, frameFBO);
		glReadBuffer(GL_COLOR_ATTACHMENT0);

		glBlitFramebuffer(0, 0, frameX, frameY, 0, 0, frameX, frameY, GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT, GL_NEAREST);

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
	printf("closing\n");
	
	glfwTerminate();
}


