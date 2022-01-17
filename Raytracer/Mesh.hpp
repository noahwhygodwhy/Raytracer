#pragma once
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
    vec3 position;
    vec3 normal;
    vec3 tangent;
    vec3 bitangent;
    vec2 texCoords;
};


struct Texture {
    vec3* data;
    string type;
    string path;
};

class Mesh
{
public:
    //vector<Mesh> children;
    vector<Vertex> vertices;
    vector<unsigned int> indices;
    mat4 transform;



    Mesh(vector<Vertex> vertices, vector<unsigned int> indices, const mat4& meshTx);
    Mesh();
    ~Mesh();
    void draw(const Shader& shader, const mat4& parentTx) const;
    //pair<vec3, vec3> getBoundingBox() const;
private:
    unsigned int VAO, VBO, EBO;

    void setupMesh();

};


#endif
