#ifndef TEXTURE_H
#define TEXTURE_H

#include <unordered_map>
#include <string>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;
using namespace std;

class Texture;

static unordered_map<string, Texture*> loadedTextures;

class Texture
{
    unsigned char* data;
    int width;
    int height;
    int nrComponents;
public:
	Texture(unsigned char* data, int width, int height, int nrComponents);
	~Texture();
    dvec3 vec3Sample(dvec2 uv);
    double doubleSample(dvec2 uv);
    std::function<double (dvec2)> getDoubleSampler() { return std::bind(&Texture::doubleSample, this, std::placeholders::_1); }
    std::function<vec3 (dvec2)> getVec3Sampler() { return std::bind(&Texture::vec3Sample, this, std::placeholders::_1); }

private:

};
Texture* loadTexture(string filename, string directory);

#endif