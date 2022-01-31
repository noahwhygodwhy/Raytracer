#ifndef TRIANGLE_H
#define TRIANGLE_H




extern uint32_t frameX;
extern uint32_t frameY;



float antiNDC(float i) {
	return (i + 1.0f) / 2.0f;
}



float triArea(vec3 a, vec3 b, vec3 c) {
	vec3 d = b - a;
	vec3 e = c - a;
	vec3 f = glm::cross(d, e);
	return length(f) / 2.0f;
}

class Triangle {
	Vertex vertices[3];
	uvec2 maxFrags;
	uvec2 minFrags;

	//uvec2 zLimits;
	bool onscreen = true;
	fvec3 U;
	fvec3 V; 
	float planeD;

public:
	vec3 normal;
	Mesh* mesh;
	Triangle() {}

	Triangle(Vertex a, Vertex b, Vertex c, Mesh* m) {
		this->mesh = m;
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;


		U = b.position - a.position;
		V = c.position - a.position;

		this->normal.x = (U.y * V.z) - (U.z * V.y);
		this->normal.y = (U.z * V.x) - (U.x * V.z);
		this->normal.z = (U.x * V.y) - (U.y * V.x);
		this->normal = glm::normalize(this->normal);
		this->planeD = glm::dot(this->normal, this->vertices[0].position);
	}
	Vertex operator[](uint32_t index) const {
		return vertices[index];
	}
private:
};
float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}
//sign() and pointInTriangle() based off of the code in this stackoverflow post https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float sign(vec2 a, vec2 b, vec2 c) {
	return ((a.x - c.x) * (b.y - c.y)) - ((b.x - c.x) * (a.y - c.y));
}
bool rayHit(const vec3& rayOrigin, const vec3& rayDirection, const Triangle& tri, vec3& barycentric) {
	
	vec3 e1 = tri[1].position - tri[0].position;
	vec3 e2 = tri[2].position - tri[0].position;


	vec3 P = rayDirection, 
	
	float denom = 1.0f/(glm::dot())


}
//bool rayHit(const vec3& rayOrigin, const vec3& rayDirection, const Triangle& tri, vec3& barycentric) {
//	printf("rayNormal: %s\n", glm::to_string(tri.normal).c_str());
//	float denom = glm::dot(tri.normal, rayDirection);
//	if (denom == 0) {
//		return false;
//	}
//
//
//	float D = glm::dot(tri.normal, tri[0].position);
//	float t = -(glm::dot(tri.normal, rayOrigin) + D) / denom;
//	/*if (t < 0) {
//		return false;
//	}*/
//	vec3 P = rayOrigin + (t * rayDirection);//P is the point in the triangle where it hits
//
//	float totalArea = triArea(tri[0].position, tri[1].position, tri[2].position);
//
//	float u = triArea(tri[2].position, tri[0].position, P) / totalArea;
//	float v = triArea(tri[0].position, tri[1].position, P) / totalArea;
//	float w = 1.0 - (v + u);// triArea(tri[1].position, tri[2].position, P) / totalArea;
//
//	//printf("uvw: %f, %f, %f\n", u, v, w);
//	if (u < 0.0f || v < 0.0f || w < 0) {
//		//printf("uvw: %f, %f, %f\n", u, v, w);
//		return false;
//	}
//	barycentric = vec3(u, v, w);
//	return true;
//
//
//	//https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
//	//needs to check if it's in the triangle, and somehow return that point, because we'll need to use it for uv and stuff
//}

bool pointInTriangle(vec2 point, Triangle tri) {
	float signab, signbc, signca;
	bool hasNeg, hasPos;

	signab = sign(point, tri[0].position, tri[1].position);
	signbc = sign(point, tri[1].position, tri[2].position);
	signca = sign(point, tri[2].position, tri[0].position);

	hasNeg = (signab < 0) || (signbc < 0) || (signca < 0);
	hasPos = (signab > 0) || (signbc > 0) || (signca > 0);

	return !(hasNeg && hasPos);
}


#endif