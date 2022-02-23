#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <string>
#include "Texture.hpp"

using namespace std;
using namespace glm;


struct materialStats {

	dvec3 color = dvec3(1.0, 0.0, 1.0);
	double ns = 100.0;//specular highlights exponent							
	double ni = 1.0; //index of refraction									
	double transparency = 0.0;//transparency
	double metalness = 0.7; //how shiny
	double smoothness = 0.5;
};

//unordered_map<string, materialStats> materials;


//TODO: finish writing up the rest of the functions
//write up the constructor that accepts only optional parameters lol
//like, bug material, and a bunch of nulls for the other functions

class Material {

public:
	Material();
	Material(
		dvec3 color,
		double ns = 100.0,
		double ni = 1.0,
		double transparency = 0.0,
		double metalness = 0.5,
		double smoothness = 0.5,
		function<dvec3(dvec2)> colorFn = nullptr,
		function<dvec3(dvec2)> normalFn = nullptr,
		function<double(dvec2)> nsFn = nullptr,
		function<double(dvec2)> niFn = nullptr,
		function<double(dvec2)> transparencyFn = nullptr,
		function<double(dvec2)> metalnessFn = nullptr,
		function<double(dvec2)> smoothnessFn = nullptr
	);
	Material(
		string premade,
		function<dvec3(dvec2)> colorFn = nullptr,
		function<dvec3(dvec2)> normalFn = nullptr,
		function<double(dvec2)> nsFn = nullptr,
		function<double(dvec2)> niFn = nullptr,
		function<double(dvec2)> transparencyFn = nullptr,
		function<double(dvec2)> metalnessFn = nullptr,
		function<double(dvec2)> smoothnessFn = nullptr
	);



	dvec3 getColor(dvec2 uv);
	dvec3 getNormal(dvec2 uv, dvec3 defaultNormal);
	double getNS(dvec2 uv);
	double getNI(dvec2 uv);
	double getTransparency(dvec2 uv);
	double getMetalness(dvec2 uv);
	double getSmoothness(dvec2 uv);
private:

	function<dvec3(dvec2)> colorFn;
	function<dvec3(dvec2)> normalFn;
	function<double(dvec2)> nsFn;
	function<double(dvec2)> niFn;
	function<double(dvec2)> transparencyFn;
	function<double(dvec2)> metalnessFn;
	function<double(dvec2)> smoothnessFn;

	dvec3 color;
	double ns;
	double ni;
	double transparency;
	double metalness;
	double smoothness;
};






#endif