#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>


using namespace std;
using namespace glm;

/*
struct Material {																//For the strauss model
	dvec3 ka = dvec3(1.0, 1.0, 1.0);											//unused
	dvec3 kd = dvec3(1.0, 1.0, 1.0);											//chromatic
	dvec3 ks = dvec3(1.0, 1.0, 1.0);											//unused
	double ns = 100.0;//specular highlights exponent							//specular highlight
	double ni = 1.0; //index of refraction										//ior
	double opaqueness = 1.0; //transparency, 0 is transparent, 1 is opaque		//transparency
	double metalness = 0.7; //how shiny											//metalness
};

unordered_map<string, Material> materials = {




	{"Brass", Material(dvec3(0.329412,0.223529,0.027451), dvec3(0.780392,0.568627,0.113725), dvec3(0.992157,0.941176,0.807843), 27.8974, 1.0, 1.0, 0.9)},
	{"Bronze", Material(dvec3(0.2125,0.1275,0.054), dvec3(0.714,0.4284,0.18144), dvec3(0.393548,0.271906,0.166721), 25.6, 1.0, 1.0, 0.9)},
	{"Polished Bronze", Material(dvec3(0.25,0.148,0.06475), dvec3(0.4,0.2368,0.1036), dvec3(0.774597,0.458561,0.200621), 76.8, 1.0, 1.0, 1.0)},
	{"Chrome", Material(dvec3(0.25,0.25,0.25), dvec3(0.4,0.4,0.4), dvec3(0.774597,0.774597,0.774597), 76.8, 1.0, 1.0)},
	{"Copper", Material(dvec3(0.19125,0.0735,0.0225), dvec3(0.7038,0.27048,0.0828), dvec3(0.256777,0.137622,0.086014), 12.8, 1.0, 1.0, 0.9)},
	{"Polished Copper", Material(dvec3(0.2295,0.08825,0.0275), dvec3(0.5508,0.2118,0.066), dvec3(0.580594,0.223257,0.0695701), 51.2, 1.0, 1.0, 1.0)},
	{"Gold", Material(dvec3(0.24725,0.1995,0.0745), dvec3(0.75164,0.60648,0.22648), dvec3(0.628281,0.555802,0.366065), 51.2, 1.0, 1.0, 0.9)},
	{"Polished Gold", Material(dvec3(0.24725,0.2245,0.0645), dvec3(0.34615,0.3143,0.0903), dvec3(0.797357,0.723991,0.208006), 83.2, 1.0, 1.0, 1.0)},
 {"Pewter", Material(dvec3(0.105882, 0.058824, 0.113725), dvec3(0.427451, 0.470588, 0.541176), dvec3(0.333333, 0.333333, 0.521569), 9.84615, 1.0, 1.0, 0.5)},
{ "Silver", Material(dvec3(0.19225,0.19225,0.19225), dvec3(0.50754,0.50754,0.50754), dvec3(0.508273,0.508273,0.508273), 51.2, 1.0, 1.0, 0.9) },
{ "Polished Silver", Material(dvec3(0.23125,0.23125,0.23125), dvec3(0.2775,0.2775,0.2775), dvec3(0.773911,0.773911,0.773911), 89.6, 1.0, 1.0, 1.0) },
{ "Emerald", Material(dvec3(0.0215,0.1745,0.0215), dvec3(0.07568,0.61424,0.07568), dvec3(0.633,0.727811,0.633), 76.8, 1.0, 0.55, 0.5) },
{ "Jade", Material(dvec3(0.135,0.2225,0.1575), dvec3(0.54,0.89,0.63), dvec3(0.316228,0.316228,0.316228), 12.8, 1.0, 0.95) },
{ "Obsidian", Material(dvec3(0.05375,0.05,0.06625), dvec3(0.18275,0.17,0.22525), dvec3(0.332741,0.328634,0.346435), 38.4, 1.0, 0.82) },
{ "Pearl", Material(dvec3(0.25,0.20725,0.20725), dvec3(1.0,0.829,0.829), dvec3(0.296648,0.296648,0.296648), 11.264, 1.0, 0.922) },
{ "Ruby", Material(dvec3(0.1745,0.01175,0.01175), dvec3(0.61424,0.04136,0.04136), dvec3(0.727811,0.626959,0.626959), 76.8, 1.0, 0.55) },
{ "Turquoise", Material(dvec3(0.1,0.18725,0.1745), dvec3(0.396,0.74151,0.69102), dvec3(0.297254,0.30829,0.306678), 12.8, 1.0, 0.8) },
{ "Black Plastic", Material(dvec3(0.0,0.0,0.0), dvec3(0.01,0.01,0.01), dvec3(0.50,0.50,0.50), 32, 1.0, 1.0, 0.0) },
{ "Black Rubber", Material(dvec3(0.02,0.02,0.02), dvec3(0.01,0.01,0.01), dvec3(0.4,0.4,0.4), 10, 1.0, 1.0, 0.0) },
{ "Bug", Material(dvec3(1.0, 0.0, 1.0), dvec3(1.0, 0.0, 1.0), dvec3(1.0, 0.0, 1.0), 1, 1.0, 1.0, 0.0) },
{ "Glass", Material(dvec3(1.0, 1.0, 1.0), dvec3(1.0, 1.0, 1.0), dvec3(1.0, 1.0, 1.0), 10, 1.52, 0.0, 0.5) },
};
*/

struct Material {																
	dvec3 kd = dvec3(1.0, 1.0, 1.0);											
	double ns = 100.0;//specular highlights exponent							
	double ni = 1.0; //index of refraction									
	double transparency  = 0.0;//transparency
	double metalness = 0.7; //how shiny
	double smoothness = 0.5;
};
unordered_map<string, Material> materials = {
	{"Glass", Material(dvec3(1.0, 1.0, 1.0), 50.0, 1.54, 0.95, 0.5, 1.0)},
	{"Bug", Material(dvec3(1.0, 0.0, 1.0), 1.0, 1.0, 0.0, 0.0, 0.0)},
	{"Copper", Material(dvec3(0.7038,0.27048,0.0828), 50.0, 1.0, 0.0, 0.5, 0.5)}

};

#endif