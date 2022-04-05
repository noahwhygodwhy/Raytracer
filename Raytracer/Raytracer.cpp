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

#define MAX_SHAPES 50
#define MAX_MATERIALS 10

#define MAX_PATH 500
//#define OUTPUTPASSES 50
#define OUTPUTFRAMES 189
//#define EVERYFRAME INFINITY
#define CONCURRENT_FOR
#define KDTRACE
//#define CIN
//#define PPCIN

#define PIXEL_MULTISAMPLE_N 1
#define MONTE_CARLO_SAMPLES 1000


//#define BASIC_BITCH


bool prd = false; //print debuging for refraction

uint32_t frameX = 1000;
uint32_t frameY = 1000;
double frameRatio = double(frameX) / double(frameY);


//dvec3 clearColor(0.21, 0.78, 0.95);
fvec4 clearColor(0.0, 0.0, 0.0, 1.0);

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


/*HitResult shootRay(const Ray& ray, const FrameInfo& fi) {

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
}*/




/*dvec3 lightingFunction(const Ray& originalRay, const Ray& lightRay, const HitResult& minRayResult, const double attenuation, const Material& mat, const dvec3& lightColor) {

	dvec3 lightReflectVector = glm::normalize((glm::dot(lightRay.direction, minRayResult.normal) * 2.0 * minRayResult.normal) - lightRay.direction);
	dvec3 H = glm::normalize(lightRay.direction + originalRay.inverseDirection);

	double spec = glm::pow(glm::max(0.0, glm::dot(lightReflectVector, originalRay.inverseDirection)), mat.getNS(minRayResult.uv));//to the specular exponent
	dvec3 specular = (lightColor * spec) / attenuation;
	double diff = glm::max(0.0, glm::dot(minRayResult.normal, lightRay.direction));
	dvec3 diffuse = (mat.getColor(minRayResult.uv) * lightColor * diff) / attenuation;
	return (diffuse + specular);
}*/

/*dvec3 randomHemisphereVector(double u1, double u2)
{
	double r = glm::sqrt(1.0f - u1 * u1);
	double phi = 2 * glm::pi<double>() * u2;

	return dvec3(cos(phi) * r, sin(phi) * r, u1);
}*/

/*dvec3 ggx(const Ray& incomingRay,
	const Ray& outgoingRay,
	const HitResult&
	minRayResult,
	dvec3 downstreamRadiance,
	Material* mat) {


}*/

/*dvec3 pathTrace(const Ray& ray, const FrameInfo& fi, double currentIOR = 1.0, uint32_t layer = 0) {
	if (prd)printf("\n======================\npath tracing layer %u\n", layer);

	if (layer > MAX_PATH) {
		return dvec3(clearColor) * 0.1;
	}
	HitResult minRayResult = shootRay(ray, fi);

	if (minRayResult.shape == NULL) {
		return dvec3(clearColor) * 0.1;
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
		else {

			if (reflectanceDecider < smooth) {
				if (prd)printf("case2\n");
				newRayDir = -rotate(ray.direction, glm::pi<double>(), (dvec3)minRayResult.normal);
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
}*/






/*void addModel(vector<Shape*>& shapes, string modelName, dvec3 pos = dvec3(0.0), dvec3 rot = dvec3(0.0, 0.0, 0.0)) {
	Model m(modelName, pos, rot);
	shapes.insert(shapes.end(), m.children.begin(), m.children.end());
}*/




AABB redoAABBs(vector<Sphere>& shapes) {

	AABB toReturn = { fvec4(0.0), fvec4(0.0) };
	for (Sphere& s : shapes) {
		s.shape.boundingBox.min = s.origin - fvec4(s.radius, s.radius, s.radius, 0.0f);
		s.shape.boundingBox.max = s.origin + fvec4(s.radius, s.radius, s.radius, 0.0f);
		toReturn.min = glm::min(toReturn.min, s.shape.boundingBox.min);
		toReturn.max = glm::max(toReturn.max, s.shape.boundingBox.max);
	}

	return toReturn;
}



int main()
{


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
	int j = 0;
	printf("\tDevice #%d:\n", j);
	//size_t size;
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

	clGetDeviceInfo(devices[j], CL_DEVICE_ADDRESS_BITS, sizeof(ui), &ui, NULL);
	printf("\t\tDevice Address Bits = %d\n", ui);

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
		
	

	cl_device_id device = devices[0];
	cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &status);
	//CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
	cl_command_queue_properties* qProperties =new cl_command_queue_properties();

	cl_command_queue cmdQueue = clCreateCommandQueueWithProperties(context, device, qProperties, &status);

	cl_mem clOtherData = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(OtherData), NULL, &status);
	cl_mem clShapes = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Shape) * MAX_SHAPES, NULL, &status);
	cl_mem clMaterials = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Material) * MAX_MATERIALS, NULL, &status);

	cl_mem clFrameBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, frameX * frameY * sizeof(frameBuffer[0]), NULL, &status);

	cl_mem clRandomBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, frameX * frameY * sizeof(uint64_t), NULL, &status);

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
	printf("build program %i\n", status);
	if (status != CL_SUCCESS)
	{ // retrieve and print the error messages:
		size_t size;
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size);
		cl_char* log = new cl_char[size];
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, size, log, NULL);
		printf("clBuildProgram failed:\n%s\n", log);
		exit(-1);
	}

	cl_kernel kernel = clCreateKernel(program, "render", &status );
	printf("made kernal %i\n", status);

	status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &clOtherData);
	printf("set args 0 %i\n", status);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &clShapes);
	printf("set args 1 %i\n", status);
	status = clSetKernelArg(kernel, 2, sizeof(cl_mem), &clMaterials);
	printf("set args 2 %i\n", status);
	status = clSetKernelArg(kernel, 3, sizeof(cl_mem), &clFrameBuffer);
	printf("set args 3 %i\n", status);
	status = clSetKernelArg(kernel, 4, sizeof(cl_mem), &clRandomBuffer);
	printf("set args 4 %i\n", status);



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

	frameBuffer = new fvec4[frameX * frameY]();
	drawBuffer = new fvec4[frameX * frameY]();


	mt19937_64 numGen;
	randomBuffer = new uint64_t[frameX * frameY]();
	for (size_t i = 0; i < frameX * frameY; i++) {
		randomBuffer[i] = numGen();
		//printf("random buffer %i: %zu\n", i, randomBuffer[i]);
	}
	//exit(0);
	

	//initialize textured that the framebuffer gets written to display it on a triangle
	unsigned int frameTexture;
	glGenTextures(1, &frameTexture);
	glBindTexture(GL_TEXTURE_2D, frameTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);



	



	//TODO: two spheres makes it hard. NEED to do the octtree


	//Material checkers("PlainWhiteTees");
	//Material* checkers = new Material("PlainWhiteTees", dvec3(0.0, 0.0, 0.0), ryCheckers10x10);

	//Material* white = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 1.0, 1.0));
	//Material* red = new Material("Red");
	//Material* green = (new Material("PlainWhiteTees"))->setColor(dvec3(0.0, 1.0, 0.0));
	//Material* blue = (new Material("PlainWhiteTees"))->setColor(dvec3(0.0, 0.0, 1.0));
	//Material* yellow = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 1.0, 0.0));
	//Material* purple = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 0.0, 1.0));
	//Material* teal = (new Material("PlainWhiteTees"))->setColor(dvec3(1.0, 0.0, 1.0));
	//Material* glass = new Material("Glass");
	//Material* mirrorA = new Material("Mirror");
	//Material* mirrorB = new Material("MirrorB");



	vector<Material> materials;
	materials.reserve(MAX_MATERIALS);
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f));
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(1.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f));
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 0.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f));
	materials.push_back(Material(fvec4(1.0f, 1.0f, 1.0f, 0.0f), fvec4(0.0f, 0.0f, 1.0f, 0.0f), 10.0f, 1.0f, 0.0f, 0.0f, 0.0f));


	vector<Sphere> shapes;
	shapes.reserve(MAX_SHAPES);
	shapes.push_back(Sphere(fvec4(0.0f, 3.0f, 0.0f, 0.0f), Shape(AABB(), 0u), 3.0f));
	shapes.push_back(Sphere(fvec4(0.0f, -1000.0f, 0.0f, 0.0f), Shape(AABB(), 2u), 1000.0f));
	shapes.push_back(Sphere(fvec4(6.0f, 6.0f, 6.0f, 0.0f), Shape(AABB(), 1u), 1.0f));
	shapes.push_back(Sphere(fvec4(-6.0f, 6.0f, 6.0f, 0.0f), Shape(AABB(), 3u), 1.0f));



	uint32_t frameCounter = 0;
	float frameTimes[30](0);
	int lastSecondFrameCount = -1;

	uint32_t fps = 30;

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
		//printf("that frame took %f seconds\n", deltaTime);
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




		//printf("scene bounding min: %s, max:%s\n", glm::to_string(sceneBounding.min).c_str(), glm::to_string(sceneBounding.max).c_str());
		//fi.kdTree = buildKDTree(fi.shapes, sceneBounding);

		constexpr double mypi = glm::pi<double>();
		clearBuffers();


		fvec3 eye = fvec3(sin(currentFrame*3.0f) * 15, 7, cos(currentFrame*3.0f) * 15);

		//fvec3 eye = fvec3(0.0f, 7.0f, 15.0f);

		fvec3 lookat = fvec3(0.0, 3.0, 0.0);
		//lookat = vec3(0.0, -5.0, 15);

		//printf("looking at %s\n", glm::to_string(lookat).c_str());
		fvec3 camForward = glm::normalize(lookat - eye);
		fvec3 camUp = glm::normalize(fvec3(0.0, 1, 0.0));
		fvec3 camRight = glm::cross(camForward, camUp);
		camUp = glm::cross(camRight, camForward);

		//fi.camPosition = eye;


		float viewPortHeight = 2.0f;
		float viewPortWidth = viewPortHeight * frameRatio;

		float fov = 90;
		float focal = (viewPortHeight / 2.0) / glm::tan(radians(fov / 2.0));
		//qua rotQuat = glm::rotation(dvec3(0.0, 0.0, -1.0), camForward);
		
		AABB sceneBounding = redoAABBs(shapes);



		Sphere asdfdsa = shapes.at(0);

		RAND_MAX;


		auto randSeed = numGen();
		

		OtherData otherData = { clearColor, fvec4(eye, 0.0f), fvec4(camRight, 0.0f), fvec4(camUp, 0.0f), fvec4(camForward, 0.0), randSeed, focal, currentFrame, 100u, uint(shapes.size()), MONTE_CARLO_SAMPLES };



		cl_event* waitAfterWrites = new cl_event[4];

		status = clEnqueueWriteBuffer(cmdQueue, clOtherData, CL_FALSE, 0, sizeof(OtherData), &otherData, 0, NULL, &waitAfterWrites[0]);
		//printf("enqueue clOtherData %i\n", status);
		status = clEnqueueWriteBuffer(cmdQueue, clShapes, CL_FALSE, 0, sizeof(Shape) * MAX_SHAPES, shapes.data(), 0, NULL, &waitAfterWrites[1]);
		//printf("enqueue clShapes %i\n", status);
		status = clEnqueueWriteBuffer(cmdQueue, clMaterials, CL_FALSE, 0, sizeof(Material) * MAX_MATERIALS, materials.data(), 0, NULL, &waitAfterWrites[2]);
		//printf("enqueue clMaterials %i\n", status);
		status = clEnqueueWriteBuffer(cmdQueue, clRandomBuffer, CL_FALSE, 0, sizeof(uint64_t) * frameX*frameY, randomBuffer, 0, NULL, &waitAfterWrites[3]);
		//printf("enqueue clMaterials %i\n", status);

		//printf("sizeof sphere on host: %zu\n", sizeof(Sphere));
		//printf("offset on host: %zu\n", offsetof(Sphere, radius));

		

		//printf("origin:");
		//printf("%s\n", glm::to_string(x.origin).c_str());
		//printf("%s\n", glm::to_string(x.shape.boundingBox.min).c_str());
		//printf("%s\n", glm::to_string(x.shape.boundingBox.max).c_str());
		//printf("matidx: %u\n", x.shape.matIdx);

		

		cl_event* waitAfterProcessing = new cl_event;


		status = clEnqueueNDRangeKernel(cmdQueue, kernel, 2, NULL, globalWorkSize, NULL, 4, waitAfterWrites, &waitAfterProcessing[i]);
		//for (int i = 0; i < MONTE_CARLO_SAMPLES; i++) {
		//	//printf("enqueue rangekernal %i\n", status);
		//}



		status = clEnqueueReadBuffer(cmdQueue, clFrameBuffer, CL_TRUE, 0, frameX * frameY * sizeof(frameBuffer[0]), frameBuffer, 1, waitAfterProcessing, NULL);
		//printf("enqueue read %i\n", status);
		


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, frameX, frameY, 0, GL_RGBA, GL_FLOAT, frameBuffer);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDrawArrays(GL_TRIANGLES, 0, 3);
		glfwSwapBuffers(window);
		processInput(window);
		glfwPollEvents();



		/*
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
		}*/




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


