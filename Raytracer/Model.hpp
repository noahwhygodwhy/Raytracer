#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <string>
#include <map>

#include <glm/glm.hpp>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "Mesh.hpp"
#include "Shader.hpp"

using namespace std;
using namespace glm;

//class Model;





class Model
{

    static map<string, Model*> loadedModels;
    unsigned int boundingVBO;


public:

    Model(string modelname);
    Model();
    ~Model();
    mat4 transform;
    vector<Mesh> children;
private:

    mat4 originalTransform;
    string directory;
    pair<vec3, vec3> boundingBox;
    void loadModel(string path);
    vector<Mesh> processNode(aiNode* node, const aiScene* scene, const mat4& parentTx, float scaleFactor);
    Mesh processMesh(aiMesh* mesh, const aiScene* scene, const mat4& nodeTx);
    vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName);

};

#endif
