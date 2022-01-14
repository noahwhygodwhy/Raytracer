#include <iostream>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace glm;



const int screenX = 800;
const int screenY = 600;

vec4 clearColor(0.0f, 0.0f, 0.5f, INFINITY);

struct Vertex {
	vec3 position;
	vec3 normal;
	vec3 color;
};

struct Triangle {
	Vertex a;
	Vertex b;
	Vertex c;
	Vertex& operator[](int index) {
		switch (index) {
			case 0:
				return a;
				break;
			case 1:
				return b;
				break;
			case 2:
				return c;
				break;
			default:
				cout << "Array index out of bound, exiting";
				exit(0);
		}
	}
};



//TODO: do i need to transform the normals as well..?
//Applies the projection and view matrix to the point
Vertex transformIt(const mat4& projection, const mat4& view, Vertex v) {
	v.position = projection * view * vec4(v.position, 1.0f);
}


//Idk if this is even what clipping is. Do i want to split triangles evenf or something this basic? I think it's more an optimization
//checks if the triangle can be clipped
bool clipTest(Triangle t) {
	return true;//TODO:
}


int main()
{
	vector<Vertex> vertices = {
		{vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.2f, 0.0f)},
		{vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.2f, 0.0f)},
		{vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.2f, 0.0f)},
	};

	float ratio = float(screenX) / float(screenY);

	mat4 projection = glm::perspective(radians(70.0f), ratio, 0.1f, 100.0f);
	mat4 view = glm::lookAt(vec3(0.0f, 0.0f, 5.0f), vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));

	vector<Triangle> triangles;


	for (int i = 0; i < vertices.size(); i += 3) {
		triangles.push_back({
				transformIt(projection, view, vertices[i + 0]),
				transformIt(projection, view, vertices[i + 1]),
				transformIt(projection, view, vertices[i + 2])
			});
	}
	

	vec4* frameBuffer = new vec4[screenX * screenY]();
	std::fill_n(frameBuffer, screenX * screenY, clearColor);

	for (const Triangle& tri : triangles) {
		//Rasterize
	}


}