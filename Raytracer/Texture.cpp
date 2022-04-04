//#include "Texture.hpp"
//
//
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//
//
//
//Texture* loadTexture(string filename, string directory) {
//    if (loadedTextures.contains(filename)) {
//        return loadedTextures.at(filename);
//    }
//    string fullFilePath = directory + '/' + filename;
//    unsigned int textureID;
//
//    int width, height, nrComponents;
//    unsigned char* data = stbi_load(fullFilePath.c_str(), &width, &height, &nrComponents, 3); //requires 3 components
//
//    if (!data) {
//        printf("error loading texture %s\n", fullFilePath.c_str());
//        exit(-1);
//    }
//
//    loadedTextures.insert({ filename, new Texture(data, width, height, nrComponents) });
//    return loadedTextures.at(filename);
//
//}
//
//double ctd(unsigned char x) {
//
//    return double(x) / 255.0;
//}
//
//
//
//
//
//
//dvec3 Texture::vec3Sample(dvec2 uv) {
//   // printf("\n\n==================\n");
//   // printf("uv: %s\n", glm::to_string(uv).c_str());
//    //printf("pixels: %s\n", glm::to_string(glm::floor(uv * dvec2(this->width, this->height))).c_str());
//
//    uv = glm::mod(uv, 1.0);
//
//    uvec2 pixelCoords = uvec2(glm::floor(uv * dvec2(this->width, this->height)));
//    //printf("pixelCoords: %u, %u\n", pixelCoords.x, pixelCoords.y);
//    size_t firstIndex = size_t(pixelCoords.x + (pixelCoords.y * this->width)) * 3u;
//    unsigned char r = this->data[firstIndex + 0u];
//    unsigned char g = this->data[firstIndex + 1u];
//    unsigned char b = this->data[firstIndex + 2u];
//
//
//    return dvec3(ctd(r), ctd(g), ctd(b));
//}
//double Texture::doubleSample(dvec2 uv) {
//    uv = glm::mod(uv, 1.0);
//    uvec2 pixelCoords = uvec2(glm::floor(uv * dvec2(this->width, this->height)));
//    size_t firstIndex = size_t(pixelCoords.x + (pixelCoords.y * width)) * 3u;
//    unsigned char r = this->data[firstIndex + 0u];
//    unsigned char g = this->data[firstIndex + 1u];
//    unsigned char b = this->data[firstIndex + 2u];
//    return (ctd(r) + ctd(g) + ctd(b)) / 3.0;
//}
//
//
//
//
//Texture::Texture(unsigned char* data, int width, int height, int nrComponents)
//{
//    this->data = data;
//    this->width = width;
//    this->height = height;
//    this->nrComponents = nrComponents;
//}
//
//Texture::~Texture()
//{
//}