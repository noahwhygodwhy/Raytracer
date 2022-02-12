#ifndef MESH_H
#define MESH_H

#include <glm/glm.hpp>
#include "glad/glad.h"
#include "Shader.hpp"
#include <vector>
#include <string>


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
    Vertex() {
        position = normal = tangent = bitangent = dvec3(1.0);
        texCoords = dvec2(1.0);
    }
};


struct Texture {
    uint32_t width;
    uint32_t height;
    unsigned char* data;
    string type;
    string path;

    vec3 sample(vec2 uv) const { //TODO: could this be more efficient? This gets done so many times
        
        uint32_t x = uint32_t(floor((uv.x * width) + 0.5f));
        uint32_t y = uint32_t(floor((uv.y * height) + 0.5f));
        if ((x + (y * width)) * 3 > width * height * 3) {
            printf("given uv: %f, %f\n", uv.x, uv.y);
            printf("x and y: %i, %i\n", x, y);
            printf("AHHHHHH@22222222\n");
            //printf("x: %i, %y")
        }
        uint64_t i = glm::min((x + (y * width)) * 3, width * height * 3);
        
        return vec3(float(data[i + 0]) / 255.0f, float(data[i + 1]) / 255.0f, float(data[i + 2]) / 255.0f);
    }
};

class Mesh
{
public:
    //vector<Mesh> children;
    Texture diffuse;
    Texture specular;
    vector<Vertex> vertices;
    vector<uint32_t> indices;
    mat4 transform;
    mat4 getModelMat() const;
    Mesh(vector<Vertex> vertices, vector<uint32_t> indices, const mat4& meshTx, const Texture& diffuse, const Texture& specular);
    Mesh();
    ~Mesh();

private:
};


#endif
