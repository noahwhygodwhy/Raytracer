#include <queue>
#include <mutex>
#include <thread>

#include <condition_variable>

#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>

#include "Model.hpp"
#include "Triangle.hpp"
using namespace std;
using namespace glm;

extern uint32_t frameX;
extern uint32_t frameY;
extern fvec2 pixelOffset;

extern vec3* frameBuffer;
extern float* depthBuffer;
mutex bufferMutex;


struct fragRef {
	Triangle tri;
	uint32_t fragx;
	uint32_t fragy;

	fragRef(Triangle t, uint32_t x, uint32_t y) { 
		tri = t;
		fragx = x;
		fragy = y;
	};
};

class WorkerPool
{
	queue <fragRef*> fragQueue;
	mutex fragQueueMutex;
	mutex emptyMutex;

public:
	WorkerPool();
	~WorkerPool();
	void runLoop(int id);
	//void addJob(fragRef* ref);
	void gibQue(queue <fragRef*>& que);
	void waitForEmptyQueue();

	void stop();


private:
	vector<thread> pool;
	condition_variable jobcondition;
	condition_variable emptycondition;

};


void WorkerPool::gibQue(queue <fragRef*>& que) {
	printf("gibQue");
	this->fragQueue = que;
	jobcondition.notify_all();

}

/*void WorkerPool::addJob(fragRef* ref) {
	{
		unique_lock<mutex> lock(fragQueueMutex);
		fragQueue.push(ref);
		//printf("pushing onto %i\n", fragQueue.size());
	}
	jobcondition.notify_one();
}*/


void renderFrag(fragRef* ref) {
	fvec2 fragCoord = fvec2((floatDiv(ref->fragx, frameX) * 2.0f) - 1.0f, (floatDiv(ref->fragy, frameY) * 2.0f) - 1.0f) + pixelOffset;
	if (pointInTriangle(fragCoord, ref->tri)) {
		vec3 baryCoords = ref->tri.getBaryCoords(fragCoord);
		if (baryCoords.x >= 0 && baryCoords.y >= 0 && baryCoords.z >= 0) {
			float fragDepth = glm::dot(vec3(ref->tri[0].position.z, ref->tri[1].position.z, ref->tri[2].position.z), baryCoords);

			vec2 fragUV = mat3x2(ref->tri[0].texCoords, ref->tri[1].texCoords, ref->tri[2].texCoords) * baryCoords;
			vec3 fragColor = ref->tri.mesh->diffuse.sample(fragUV);
			{//Section that locks the buffer. Try to do stuff outside of these brackets if possible
				//std::unique_lock<std::mutex> lock(bufferMutex);
				if (fragDepth < depthBuffer[ref->fragx + (ref->fragy * frameX)]) {
					depthBuffer[ref->fragx + (ref->fragy * frameX)] = fragDepth;
					frameBuffer[ref->fragx + (ref->fragy * frameX)] = fragColor;
				}
			}
		}
	}

}

bool terminatePool = false;
void WorkerPool::runLoop(int id) {
	fragRef* Job;
	while (true) {

		/*Job = fragQueue.front();
		fragQueue.pop();*/

		 {
			//printf("thread %i is waiting for the queue to be free\n", id);
			std::unique_lock<std::mutex> lock(fragQueueMutex);
			//printf("thread %i has the frag queue\n", id);
			jobcondition.wait(lock, [this]() {
				return !fragQueue.empty() && !terminatePool;
				});
			Job = fragQueue.front();
			fragQueue.pop();
			//printf("popping off %i\n", fragQueue.size());
			lock.unlock();
			emptycondition.notify_one();
		}
		renderFrag(Job);
	}
}


WorkerPool::WorkerPool() {
	int numThreads = 1;// std::thread::hardware_concurrency();
	//printf("num_threads: %i", numThreads);
	//num_threads = 1;
	//std::vector<std::thread> threads;
	for (int i = 0; i < numThreads; i++)
	{
		pool.push_back(std::thread(&WorkerPool::runLoop, this, i));
	}
}
WorkerPool::~WorkerPool() {

}

void WorkerPool::waitForEmptyQueue() {
	{
		std::unique_lock<std::mutex> lock(emptyMutex);
		emptycondition.wait(lock, [this] {return fragQueue.empty(); });
		emptycondition.notify_one();
		lock.unlock();

	}
}
void WorkerPool::stop() {
	terminatePool = true;
}

