#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "Ray.hpp"
#include "Shape.hpp"
#include "Material.hpp"
#include "Texture.hpp"
#include "Procedural.hpp"

extern bool prd;

using namespace std;
using namespace glm;

struct Vertex {
	dvec3 position;
	dvec3 normal;
	dvec3 tangent;
	dvec3 bitangent;
	dvec2 texCoords;
	Vertex(dvec3 a) {
		position = a;
	}
	Vertex(dvec3 a, dvec2 uv) {
		position = a;
		texCoords = uv;
	}
	Vertex(dvec3 a, dvec3 n) {
		position = a;
		normal = n;
	}
	Vertex(dvec3 a, dvec3 n, dvec2 uv) {
		position = a;
		normal = n;
		texCoords = uv;
	}
	Vertex() {
		position = normal = tangent = bitangent = dvec3(1.0);
		texCoords = dvec2(1.0);
	}
};



class Triangle : public Shape {
	Vertex vertices[3];

public:
	Triangle(Vertex a, Vertex b, Vertex c, Material* material)
		: Shape(AABB(glm::min(a.position, glm::min(b.position, c.position)),
			glm::max(a.position, glm::max(b.position, c.position))), material)
	{
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
	}

	vec3 getP(const mat3x3& transformedVerts, vec3 bary) const {
		return mat3x3(vertices[0].position, vertices[1].position, vertices[2].position) * bary;
	}
	Vertex operator[](uint32_t index) const {
		return vertices[index];
	}
	bool rayHit(const Ray& ray, HitResult& hit, double currentTime)const;
	void redoAABB(double currentTime);
private:

};






#endif