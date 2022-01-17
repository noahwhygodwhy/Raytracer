#include "Mesh.hpp"

#include "glm/gtx/string_cast.hpp"

Mesh::Mesh(vector<Vertex> vertices, vector<unsigned int> indices, const mat4& meshTx) {

    printf("size of vertices: %i\n", vertices.size());
    printf("first vertex: %s\n", glm::to_string(vertices.at(0).position).c_str());
    this->vertices = vertices;
    this->indices = indices;
    this->transform = meshTx;
    this->setupMesh();
}

Mesh::Mesh() {

}


Mesh::~Mesh() {

}
/*
void Mesh::setupMesh() {

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

    // vertex positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    // vertex normals
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    // vertex texture coords
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tangent));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, bitangent));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, texCoords));

    glBindVertexArray(0);
}*/


/*void Mesh::draw(const Shader& shader, const mat4& parentTx) const {

    //printf("drawing mesh\n");
    //shader.use();
    mat4 tx = parentTx * this->transform;
    shader.setMatFour("model", tx);
    shader.setMatFour("normalMatrix", glm::transpose(glm::inverse(tx)));

    glBindVertexArray(VAO);

    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}*/

