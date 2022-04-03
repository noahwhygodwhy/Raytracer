#include "Raytracer.hpp"
#include "CookTorrance.hpp"

using namespace std;
using namespace std::filesystem;
using namespace glm;


#define MAX_PATH 500
//#define OUTPUTPASSES 50
#define OUTPUTFRAMES 189
#define EVERYFRAME INFINITY
#define CONCURRENT_FOR
#define KDTRACE
//#define CIN
//#define PPCIN

#define PIXEL_MULTISAMPLE_N 4
#define MONTE_CARLO_SAMPLES 1000


//#define BASIC_BITCH


bool prd = false; //print debuging for refraction

uint32_t frameX = 1000;
uint32_t frameY = 1000;
double frameRatio = double(frameX) / double(frameY);


//dvec3 clearColor(0.21, 0.78, 0.95);
dvec3 clearColor(0.0, 0.0, 0.0);

double deltaTime = 0.0f;	// Time between current frame and last frame
double lastFrame = 0.0f; // Time of last frame
string saveFileDirectory = "";


constexpr double bias = 1e-4;


void frameBufferSizeCallback(GLFWwindow* window, uint64_t width, uint64_t height) {
	glViewport(0, 0, GLsizei(width), GLsizei(height));
}

void clearBuffers() {
	for (uint64_t i = 0; i < frameX * frameY; i++) {
		frameBuffer[i] = clearColor;
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


	double closeness = glm::dot(hitNormal, incidentVector);
	double prevIOR = 1.0;
	double newIOR = objectIOR;


	if (!entering) {
		hitNormal = -hitNormal;
		swap(prevIOR, newIOR);
	}

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
		return incidentVector;
		//exit(-1);
	}

	double cosA2 = glm::sqrt(1.0 - (sinA2 * sinA2));

	if (cosA1 < 0.0) {
		cosA2 *= -1.0;
	}

	return glm::normalize(trueReflectDir);
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

dvec3 ggx(const Ray& incomingRay,
	const Ray& outgoingRay,
	const HitResult&
	minRayResult,
	dvec3 downstreamRadiance,
	Material* mat) {


}

dvec3 pathTrace(const Ray& ray, const FrameInfo& fi, double currentIOR = 1.0, uint32_t layer = 0) {
	if (prd)printf("\n======================\npath tracing layer %u\n", layer);

	if (layer > MAX_PATH) {
		return clearColor * 0.1;
	}
	HitResult minRayResult = shootRay(ray, fi);

	if (minRayResult.shape == NULL) {
		return clearColor * 0.1;
	}
#ifdef BASIC_BITCH
	return minRayResult.shape->mat->getColor(minRayResult.uv);
#endif

	dvec2 uv = minRayResult.uv;
	Material* mat = minRayResult.shape->mat;

	double trans = mat->getTransparency(uv);
	double smooth = mat->getSmoothness(uv);
	double metal = mat->getMetalness(uv);

	double transparencyDecider = randDubTwo();
	double reflectanceDecider = randDubTwo();
	double specularDecider = randDubTwo();


	dvec3 newRayDirection;
	dvec3 downstreamRadiance;


	double hitAngle = glm::acos(glm::dot(minRayResult.normal, ray.inverseDirection));
	bool entering = hitAngle < (glm::pi<double>() / 2.0);
	bool internalOnly;

	dvec3 newRayDir;
	dvec3 newRayPos;

	if (prd)printf("transD: %f, reflect: %f\n", transparencyDecider, reflectanceDecider);


	dvec3 thisRadiance = mat->getEmission();
	if (thisRadiance == dvec3(0.0, 0.0, 0.0)) {

		if (transparencyDecider < trans) {
			if (prd)printf("case1\n");

			newRayDir = glm::normalize(getRefractionRay(glm::normalize(minRayResult.normal), glm::normalize(ray.direction), mat->getNI(minRayResult.uv), entering, internalOnly));
			newRayPos = minRayResult.position + (minRayResult.normal * (entering ? -1.0 : 1.0) * bias);

			Ray reflectionRay(newRayPos, newRayDir);


			if (prd)printf("radiance is 0\n");

			thisRadiance = pathTrace(reflectionRay, fi, mat->getNI(uv), layer + 1);

		}
		/*else if (reflectanceDecider < smooth) {
			if (prd)printf("case2\n");
			newRayDir = -rotate(ray.direction, glm::pi<double>(), minRayResult.normal);
			newRayPos = minRayResult.position + (minRayResult.normal * bias);
		}*/
		else {

			if (reflectanceDecider < smooth) {
				if (prd)printf("case2\n");
				newRayDir = -rotate(ray.direction, glm::pi<double>(), minRayResult.normal);
				newRayPos = minRayResult.position + (minRayResult.normal * bias);
			}
			else {
				if (prd)printf("case3\n");
				newRayDir = glm::normalize(randomHemisphericalVector(minRayResult.normal));
				newRayPos = minRayResult.position + (minRayResult.normal * bias);
			}

		

			if (prd)printf("reflectionRay: %s, %s\n", glm::to_string(newRayPos).c_str(), glm::to_string(newRayDir).c_str());

			Ray reflectionRay(newRayPos, newRayDir);


			if(prd)printf("radiance is 0\n");
		
			downstreamRadiance = pathTrace(reflectionRay, fi, mat->getNI(uv), layer + 1);

			//thisRadiance += ggx(ray, reflectionRay, minRayResult, downstreamRadiance, mat);


			dvec3 kS = dvec3(0);

			dvec3 F0a = dvec3(glm::abs((1.0 - mat->getNI(uv)) / (1.0 + mat->getNI(uv))));
			F0a = F0a * F0a;
			dvec3 matC = mat->getColor(uv);

			dvec3 F0 = glm::mix(F0a, matC, metal);

			if(prd)printf("back on layer %u\n", layer);
			if (prd)printf("F0: %s\n", glm::to_string(F0).c_str());




			//dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, minRayResult.shape->mat, currentIOR, F0, kS);
			dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, mat, currentIOR, F0, kS);




			if(prd)printf("ks: %f\n", kS);

	//blinn phong (ish)
			dvec3 R = glm::rotate(newRayDir, glm::pi<double>(), minRayResult.normal);
			dvec3 V = ray.origin;
			dvec3 N = minRayResult.normal;
			dvec3 L = newRayDir;

			dvec3 H = (L + V) / glm::length(L + V);

			double diff = glm::dot(L, N);
			double spec = glm::dot(N, H);

			//double ks = smooth;
			dvec3 kD = (1.0 - kS) * (1.0 - metal);;

			dvec3 diffuse = diff * downstreamRadiance * mat->getColor(uv);
			//dvec3 specular = kS * (pow(spec, mat->getNI(uv))) * downstreamRadiance * mat->getColor(uv);

			//thisRadiance += (downstreamRadiance* mat->getColor(uv));

			thisRadiance += (ctSpecular)+(kD * diffuse);
		}
	}

	if (prd)printf("thisradiance: %s\n", glm::to_string(thisRadiance).c_str());

	return thisRadiance;
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
	//cl_device_id device_id = NULL;
	//cl_context context = NULL;
	//cl_command_queue command_queue = NULL;
	//cl_mem memobj = NULL;
	//cl_program program = NULL;
	//cl_kernel kernel = NULL;
	//cl_platform_id platform_id = NULL;
	//cl_uint ret_num_devices;
	//cl_uint ret_num_platforms;
	//cl_int ret;

	const int NUM_ELEMENTS = 10;
	int* firstIn = new int[NUM_ELEMENTS]();
	for (int i = 0; i < NUM_ELEMENTS; i++) {
		firstIn[i] = i;
	}
	int* firstOut = new int[NUM_ELEMENTS]();

	size_t dataSize = NUM_ELEMENTS * sizeof(int);




	cl_uint numPlatforms;
	cl_int status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetPlatformIDs failed (1)\n");
	printf("Number of Platforms = %d\n", numPlatforms);
	cl_platform_id* platforms = new cl_platform_id[numPlatforms];
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetPlatformIDs failed (2)\n");
	cl_uint numDevices;
	cl_device_id* devices;



	int i = 0;
	
	printf("Platform #%d:\n", i);
	size_t size;
	char* str;
	clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &size);
	str = new char[size];
	clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, size, str, NULL);
	printf("\tName = '%s'\n", str);
	delete[] str;
	clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 0, NULL, &size);
	str = new char[size];
	clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, size, str, NULL);
	printf("\tVendor = '%s'\n", str);
	delete[] str;
	clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 0, NULL, &size);
	str = new char[size];
	clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, size, str, NULL);
	printf("\tVersion = '%s'\n", str);
	delete[] str;
	clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, 0, NULL, &size);
	str = new char[size];
	clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, size, str, NULL);
	printf("\tProfile = '%s'\n", str);
	delete[] str;
	// find out how many devices are attached to each platform and get their ids:
	status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetDeviceIDs failed (2)\n");
	devices = new cl_device_id[numDevices];
	status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
	if (status != CL_SUCCESS)
		fprintf(stderr, "clGetDeviceIDs failed (2)\n");
	for (int j = 0; j < (int)numDevices; j++)
	{
		printf("\tDevice #%d:\n", j);
		size_t size;
		cl_device_type type;
		cl_uint ui;
		size_t sizes[3] = { 0, 0, 0 };
		clGetDeviceInfo(devices[j], CL_DEVICE_TYPE, sizeof(type), &type, NULL);
		printf("\t\tType = 0x%04x = ", type);
		switch (type)
		{
		case CL_DEVICE_TYPE_CPU:
			printf("CL_DEVICE_TYPE_CPU\n");
			break;
		case CL_DEVICE_TYPE_GPU:
			printf("CL_DEVICE_TYPE_GPU\n");
			break;
		case CL_DEVICE_TYPE_ACCELERATOR:
			printf("CL_DEVICE_TYPE_ACCELERATOR\n");
			break;
		default:
			printf("Other...\n");
			break;
		}
		clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR_ID, sizeof(ui), &ui, NULL);
		printf("\t\tDevice Vendor ID = 0x%04x\n", ui);
		clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(ui), &ui, NULL);
		printf("\t\tDevice Maximum Compute Units = %d\n", ui);
		clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(ui), &ui, NULL);
		printf("\t\tDevice Maximum Work Item Dimensions = %d\n", ui);
		clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(sizes), sizes, NULL);
		printf("\t\tDevice Maximum Work Item Sizes = %d x %d x %d\n", sizes[0], sizes[1], sizes[2]);
		clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size), &size, NULL);
		printf("\t\tDevice Maximum Work Group Size = %d\n", size);
		clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(ui), &ui, NULL);
		printf("\t\tDevice Maximum Clock Frequency = %d MHz\n", ui);

		size_t extensionSize;
		clGetDeviceInfo(devices[j], CL_DEVICE_EXTENSIONS, 0, NULL, &extensionSize);
		char* extensions = new char[extensionSize];
		clGetDeviceInfo(devices[j], CL_DEVICE_EXTENSIONS, extensionSize, extensions, NULL);
		fprintf(stderr, "\nDevice Extensions:\n");
		for (int i = 0; i < (int)strlen(extensions); i++)
		{
			if (extensions[i] == ' ')
				extensions[i] = '\n';
		}
		fprintf(stderr, "%s\n", extensions);
		delete[] extensions;
		
	}

	cl_device_id device = devices[0];
	cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &status);

	cl_command_queue_properties* qProperties =new cl_command_queue_properties(CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);

	cl_command_queue cmdQueue = clCreateCommandQueueWithProperties(context, device, qProperties, &status);
	//cl_command_queue cmdQueue = clCreateCommandQueue(context, device, 0, &status);

	cl_mem dA = clCreateBuffer(context, CL_MEM_READ_ONLY, dataSize, NULL, &status);
	cl_mem dB = clCreateBuffer(context, CL_MEM_WRITE_ONLY, dataSize, NULL, &status);

	status = clEnqueueWriteBuffer(cmdQueue, dA, CL_FALSE, 0, dataSize, firstIn, 0, NULL, NULL);

	long unsigned int length;
	ifstream stream("Renderer.cl", ios::in | ios::ate | ios::binary);
	//stream.seekg(0, ios::end);
	length = long unsigned int(stream.tellg());
	stream.seekg(0, ios::beg);
	char* shaderSource = new char[length + 1];
	shaderSource[length] = '\0';
	const char** shaderSourcesArray = new const char* [1];
	stream.read(shaderSource, length);
	printf("shader source: %s\n", shaderSource);

	shaderSourcesArray[0] = shaderSource;

	cl_program program = clCreateProgramWithSource(context, 1, shaderSourcesArray, NULL, &status);

	const char* options = { "" };
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

	cl_kernel kernel = clCreateKernel(program, "square", &status );
	status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &dA);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &dB);


	size_t globalWorkSize[3] = { NUM_ELEMENTS, 1, 1 };
	size_t localWorkSize[3] = { 1, 1, 1 };

	status = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
	
	status = clEnqueueReadBuffer(cmdQueue, dB, CL_TRUE, 0, dataSize, firstOut, 0, NULL, NULL);

	for (int l = 0; l < NUM_ELEMENTS; l++) {
		printf("%i -> %i\n", firstIn[l], firstOut[l]);
	}

	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clReleaseCommandQueue(cmdQueue);
	clReleaseMemObject(dA);
	clReleaseMemObject(dB);

	delete[] firstIn;
	delete[] firstOut;


	exit(0);


	
	//ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
	//ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
	//context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

	///* Create Command Queue */
	//command_queue = clCreateCommandQueueWithProperties(context, device_id, 0, &ret);

	///* Create Memory Buffer */
	//memobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 1 * sizeof(int), NULL, &ret);

	///* Create Kernel Program from the source */
	//program = clCreateProgramWithSource(context, 1, (const char**)&shaderSource,
	//	(const size_t*)&shaderSource, &ret);

	///* Build Kernel Program */
	//ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

	///* Create OpenCL Kernel */
	//kernel = clCreateKernel(program, "hello", &ret);

	///* Set OpenCL Kernel Parameters */
	//ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&memobj);

	///* Execute OpenCL Kernel */
	//ret = clEnqueueNDRangeKernel(command_queue, kernel, 0, NULL, NULL, NULL, NULL, NULL, NULL);
	//

	///* Copy results from the memory buffer */
	//ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
	//	1 * sizeof(int), firstIn, 0, NULL, NULL);

	///* Display Result */
	//printf("%i\n", *firstOut);

	/* Finalization */
	/*ret = clFlush(command_queue);
	ret = clFinish(command_queue);
	ret = clReleaseKernel(kernel);
	ret = clReleaseProgram(program);
	ret = clReleaseMemObject(memobj);
	ret = clReleaseCommandQueue(command_queue);
	ret = clReleaseContext(context);

	free(shaderSource);*/
	exit(-1);

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
	Material* checkers = new Material("PlainWhiteTees", dvec3(0.0, 0.0, 0.0), ryCheckers10x10);

	Material* white = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 1.0, 1.0));
	Material* red = new Material("Red");
	Material* green = (new Material("PlainWhiteTees"))->setColor(dvec3(0.0, 1.0, 0.0));
	Material* blue = (new Material("PlainWhiteTees"))->setColor(dvec3(0.0, 0.0, 1.0));
	Material* yellow = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 1.0, 0.0));
	Material* purple = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 0.0, 1.0));
	Material* teal = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 0.0, 1.0));
	Material* glass = new Material("Glass");
	Material* mirrorA = new Material("Mirror");
	Material* mirrorB = new Material("MirrorB");


	//Material pyt("PlainWhiteTees");


	Vertex a(dvec3(-1000.0, 0, -1000.0), dvec3(0.0, 1.0, 0.0), dvec2(0.0, 0.0));
	Vertex b(dvec3(1000.0, 0, -1000.0), dvec3(0.0, 1.0, 0.0), dvec2(1.0, 0.0));
	Vertex c(dvec3(1000.0, 0, 1000.0), dvec3(0.0, 1.0, 0.0), dvec2(1.0, 1.0));
	Vertex d(dvec3(-1000.0, 0, 1000.0), dvec3(0.0, 1.0, 0.0), dvec2(0.0, 1.0));

	/*

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



	shapes.push_back(new Triangle(a, c, b, white));
	
	shapes.push_back(new Triangle(a, d, c, white));


	double wallRadius = 100000.0;
	double roomSize = 7.0;
	double totalRad = wallRadius + roomSize;

	//shapes.push_back(new Sphere(dvec3(0.0, totalRad, 0.0), wallRadius, white, noMovement));
	//shapes.push_back(new Sphere(dvec3(0.0, -totalRad, 0.0), wallRadius, white, noMovement));
	//shapes.push_back(new Sphere(dvec3(totalRad, 0.0, 0.0), wallRadius, red, noMovement));
	//shapes.push_back(new Sphere(dvec3(-totalRad, 0.0, 0.0), wallRadius, green, noMovement));
	//shapes.push_back(new Sphere(dvec3(0.0, 0.0, totalRad), wallRadius, wg, noMovement));
	///shapes.push_back(new Sphere(dvec3(0.0, 0.0, -totalRad), wallRadius, white, noMovement));


	shapes.push_back(new Sphere(dvec3(0.0, 2.5, 0.0), 5.0, mirrorA, noMovement));


	shapes.push_back(new Sphere(dvec3(-7.0, 1.5, 7.0), 1.5, glass, noMovement));
	shapes.push_back(new Sphere(dvec3(-7.0, 1.5, -7.0), 1.5, glass, noMovement));
	shapes.push_back(new Sphere(dvec3(7.0, 1.5, 7.0), 1.5, glass, noMovement));
	shapes.push_back(new Sphere(dvec3(7.0, 1.5, -7.0), 1.5, glass, noMovement));



	shapes.push_back(new Sphere(dvec3(0.0, 0.0, 0.0), 1, new Material("PlainWhiteTees", dvec3(1.0, 0.0, 0.0)), circle0));
	shapes.push_back(new Sphere(dvec3(0.0, 0.0, 0.0), 1, new Material("PlainWhiteTees", dvec3(0.0, 1.0, 0.0)), circle1));
	shapes.push_back(new Sphere(dvec3(0.0, 0.0, 0.0), 1, new Material("PlainWhiteTees", dvec3(0.0, 0.0, 1.0)), circle2));

	//shapes.push_back(new Sphere(dvec3(0.0, 0, 0.0), 2, red, noMovement));
	//shapes.push_back(new Sphere(dvec3(0.0, 0, -4.0), 2, mirrorA, noMovement));
	


	//shapes.push_back(new Sphere(dvec3(0.0, 3, 0.0), 3, red, circle0));
	//shapes.push_back(new Sphere(dvec3(-6.0, 3, 0.0), 3, green, noMovement));
	//shapes.push_back(new Sphere(dvec3(6.0, 3, 0.0), 3, blue, noMovement));
	//shapes.push_back(new Sphere(dvec3(2.0, 5, 5), 2, glass, noMovement));
	//shapes.push_back(new Sphere(dvec3(-6.0, 8, -8.0), 3, mirrorA, noMovement));
	//shapes.push_back(new Sphere(dvec3(6.0, 8, -8.0), 3, glass, noMovement));


    //shapes.push_back(new Sphere(dvec3(1.0, -4.2, 0.0), 2.0, Material("Glass"), noMovement));

	constexpr double mypifornow = glm::pi<double>();
	//addModel(shapes, "backpack");
	//addModel(shapes, "bunny", vec3(-0.0, 0.0, -0.0));

	//vector<Light*> lights;




	FrameInfo fi;
	fi.currentTime = 0;
	fi.shapes = shapes;
	//fi.lights = lights;
	fi.kdTree = NULL;



	//printKDTree(fi.kdTree);


	int frameCounter = -1;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	uint32_t fps = 30;


	/*time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/" + saveFileDirectory).c_str());*/
#if defined(OUTPUTFRAMES)|| defined(OUTPUTPASSES)

	time_t now;
	time(&now);
	char buf[sizeof "####-##-##-##-##-##"];
	strftime(buf, sizeof buf, "%F-%H-%M-%S", gmtime(&now));
	saveFileDirectory = string(buf);


	printf("about to create directory: %s\n", saveFileDirectory.c_str());
	filesystem::create_directory(("out/"+saveFileDirectory).c_str());
	double currentFrame = double(frameCounter) / double(fps);
	printf("frame counter: %u\n", frameCounter);

#ifdef OUTPUTFRAMES
	for (uint32_t frameCounter = 0; frameCounter < OUTPUTFRAMES;) {
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

		AABB sceneBounding = redoAABBs(fi);
		//printf("scene bounding min: %s, max:%s\n", glm::to_string(sceneBounding.min).c_str(), glm::to_string(sceneBounding.max).c_str());
		fi.kdTree = buildKDTree(fi.shapes, sceneBounding);

		constexpr double mypi = glm::pi<double>();
		clearBuffers();

		//dvec3 eye = vec3(sin(currentFrame) * 8, 2, cos(currentFrame) * 8);
		dvec3 eye = dvec3(0, 35, 50);

		dvec3 lookat = vec3(0.0, 0.0, 0.0);
		//lookat = vec3(0.0, -5.0, 15);

		printf("looking at %s\n", glm::to_string(lookat).c_str());
		dvec3 camForward = glm::normalize(lookat - eye);
		dvec3 camUp = glm::normalize(vec3(0.0, 1, 0.0));
		dvec3 camRight = glm::cross(camForward, camUp);
		camUp = glm::cross(camRight, camForward);

		fi.camPosition = eye;


		double viewPortHeight = 2.0f;
		double viewPortWidth = viewPortHeight * frameRatio;

		double fov = 90;
		double focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));
		//qua rotQuat = glm::rotation(dvec3(0.0, 0.0, -1.0), camForward);



		clearBuffers();
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

		        prd = (x == 500 && y == 250) && prd;
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
				frameBuffer[j] = ((drawBuffer[j] / float(i+1)));
				frameBuffer[j].x = cbrt(frameBuffer[j].x);
				frameBuffer[j].y = cbrt(frameBuffer[j].y);
				frameBuffer[j].z = cbrt(frameBuffer[j].z);
			}
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, frameX, frameY, 0, GL_RGB, GL_FLOAT, frameBuffer);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glDrawArrays(GL_TRIANGLES, 0, 3);
			glfwSwapBuffers(window);
			processInput(window);
			glfwPollEvents();
			
#ifdef OUTPUTPASSES
			if (i % OUTPUTPASSES == 0) {
				printf("outputting pass at i %i with name %s\n", i, (std::to_string((i / OUTPUTPASSES) + 1) + ".png").c_str());
				saveImage((std::to_string((i/OUTPUTPASSES)+1) + ".png"), window);
			}
#endif
#ifdef PPCIN
			cin.get();
#endif
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


