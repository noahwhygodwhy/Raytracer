//#ifndef MATERIAL_H
//#define MATERIAL_H
//
//#include <vector>
//#include <unordered_map>
//#include <glm/glm.hpp>
//#include <glm/gtx/string_cast.hpp>
//#include <string>
//#include "Texture.hpp"
//
//using namespace std;
//using namespace glm;
//
//
//struct materialStats {
//
//	dvec3 color = dvec3(1.0, 0.0, 1.0);
//	double ns = 100.0;//specular highlights exponent							
//	double ni = 1.0; //index of refraction									
//	double transparency = 0.0;//transparency
//	double metalness = 0.7; //how shiny
//	double smoothness = 0.5;
//};
//
//class Material {
//
//public:
//	Material();
//	Material* setColor(dvec3 newColor);
//	Material(
//		dvec3 color,
//		double ns = 100.0,
//		double ni = 1.0,
//		double transparency = 0.0,
//		double metalness = 0.5,
//		double smoothness = 0.5,
//		dvec3 emission = dvec3(0.0, 0.0, 0.0),
//		function<dvec3(dvec2)> colorFn = nullptr,
//		function<dvec3(dvec2)> normalFn = nullptr,
//		function<double(dvec2)> nsFn = nullptr,
//		function<double(dvec2)> niFn = nullptr,
//		function<double(dvec2)> transparencyFn = nullptr,
//		function<double(dvec2)> metalnessFn = nullptr,
//		function<double(dvec2)> smoothnessFn = nullptr
//	);
//	Material(
//		string premade,
//		dvec3 emmision = dvec3(0.0, 0.0, 0.0),
//		function<dvec3(dvec2)> colorFn = nullptr,
//		function<dvec3(dvec2)> normalFn = nullptr,
//		function<double(dvec2)> nsFn = nullptr,
//		function<double(dvec2)> niFn = nullptr,
//		function<double(dvec2)> transparencyFn = nullptr,
//		function<double(dvec2)> metalnessFn = nullptr,
//		function<double(dvec2)> smoothnessFn = nullptr
//	);
//
//
//
//	dvec3 getColor(dvec2 uv)const;
//	dvec3 getNormal(dvec2 uv, dvec3 defaultNormal)const;
//	double getNS(dvec2 uv)const;
//	double getNI(dvec2 uv)const;
//	double getTransparency(dvec2 uv)const;
//	double getMetalness(dvec2 uv)const;
//	double getSmoothness(dvec2 uv)const;
//	dvec3 getEmission() const;
//private:
//
//	function<dvec3(dvec2)> colorFn;
//	function<dvec3(dvec2)> normalFn;
//	function<double(dvec2)> nsFn;
//	function<double(dvec2)> niFn;
//	function<double(dvec2)> transparencyFn;
//	function<double(dvec2)> metalnessFn;
//	function<double(dvec2)> smoothnessFn;
//
//	dvec3 color;
//	double ns;
//	double ni;
//	double transparency;
//	double metalness;
//	double smoothness;
//	dvec3 emission;
//};
//
//
//
//
//
//
//#endif