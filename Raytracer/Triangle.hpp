#ifndef TRIANGLE_H
#define TRIANGLE_H




extern uint32_t frameX;
extern uint32_t frameY;



float antiNDC(float i) {
	return (i + 1.0f) / 2.0f;
}

struct Triangle {
	Vertex vertices[3];
	uvec2 maxFrags;
	uvec2 minFrags;

	Mesh* mesh;

	//uvec2 zLimits;
	bool onscreen = true;

	fvec2 v0;//for barycentry coords
	fvec2 v1;//^
	float d00;
	float d01;
	float d11;
	float denom;


	Triangle() {}

	Triangle(Vertex a, Vertex b, Vertex c, Mesh* m) {
		this->mesh = m;
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;


		v0 = b.position.xy - a.position.xy;
		v1 = c.position.xy - a.position.xy;

		d00 = glm::dot(v0, v0);
		d01 = glm::dot(v0, v1);
		d11 = glm::dot(v1, v1);

		denom = (d00 * d11) - (d01 * d01);


		//printf("frame dims: %u, %u\n", frameX, frameY);

		//printf("a: %s\n", glm::to_string(a.position).c_str());
		//printf("b: %s\n", glm::to_string(b.position).c_str());
		//printf("c: %s\n", glm::to_string(c.position).c_str());

		float maxClipX = glm::max(glm::max(a.position.x, b.position.x), c.position.x);
		float minClipX = glm::min(glm::min(a.position.x, b.position.x), c.position.x);
		float maxClipY = glm::max(glm::max(a.position.y, b.position.y), c.position.y);
		float minClipY = glm::min(glm::min(a.position.y, b.position.y), c.position.y);
		//printf("maxClip %f, %f\n", maxClipX, maxClipY);
		//printf("minClip %f, %f\n", minClipX, minClipY);

		maxClipX = antiNDC(maxClipX);
		minClipX = antiNDC(minClipX);
		maxClipY = antiNDC(maxClipY);
		minClipY = antiNDC(minClipY);

		//printf("antindcmaxClip %f, %f\n", maxClipX, maxClipY);
		//printf("antindcminClip %f, %f\n", minClipX, minClipY);

		maxFrags.x = glm::min(uint32_t(floor((maxClipX * float(frameX))+0.5f)), frameX-1);
		minFrags.x = glm::max(uint32_t(floor((minClipX * float(frameX))+0.5f)), 0u);
		maxFrags.y = glm::min(uint32_t(floor((maxClipY * float(frameY))+0.5f)), frameY-1);
		minFrags.y = glm::max(uint32_t(floor((minClipY * float(frameY))+0.5f)), 0u);

		//printf("new:\n");
		//printf("maxfrags: %s\n", glm::to_string(maxFrags).c_str());
		//printf("minfrags: %s\n", glm::to_string(minFrags).c_str());


		/*maxFrags.x = uint32_t(ceil(glm::min((glm::max(glm::max(a.position.x, b.position.x), c.position.x) + 1.0f) / 2.0f * frameX, float(frameX - 1))));
		minFrags.x = uint32_t(floor(glm::max((glm::min(glm::min(a.position.x, b.position.x), c.position.x) + 1.0f) / 2.0f * frameX, 0.0f)));
		maxFrags.y = uint32_t(ceil(glm::min((glm::max(glm::max(a.position.y, b.position.y), c.position.y) + 1.0f) / 2.0f * frameY, float(frameY - 1))));
		minFrags.y = uint32_t(floor(glm::max((glm::min(glm::min(a.position.y, b.position.y), c.position.y) + 1.0f) / 2.0f * frameY, 0.0f)));*/

		//printf("old:\n");
		//printf("maxfrags: %s\n", glm::to_string(maxFrags).c_str());
		//printf("minfrags: %s\n", glm::to_string(minFrags).c_str());

		/*if (maxFrags.y >= frameY) {
			printf("max frags %i\n", maxFrags.y);
			printf("maxClipx: %f, %f\n", maxClipX, maxClipY);
			printf("minClipx: %f, %f\n", minClipX, minClipY);
		}*/

		//exit(0);

		//zLimits.x = glm::min(glm::min(a.position.z, b.position.z), c.position.z);
		//zLimits.y = glm::max(glm::max(a.position.z, b.position.z), c.position.z);

		if (minClipX > 1.0f || minClipY > 1.0f || maxClipX < -1.0f || maxClipY < -1.0f) {
			onscreen = false;
		}
		/*else {
			printf("on screen = true");
		}*/
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
float floatDiv(uint64_t a, uint64_t b) {
	return float(a) / float(b);
}
//sign() and pointInTriangle() based off of the code in this stackoverflow post https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
float sign(vec2 a, vec2 b, vec2 c) {
	return ((a.x - c.x) * (b.y - c.y)) - ((b.x - c.x) * (a.y - c.y));
}

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