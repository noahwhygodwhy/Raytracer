#include "Material.hpp"







//lol this is bad
unordered_map<string, materialStats> materials = unordered_map<string, materialStats>{
	{"Glass", materialStats(dvec3(1.0, 1.0, 1.0), 50.0, 1.54, 0.8, 0.2, 1.0)},
	{"PlainWhiteTees", materialStats(dvec3(1.0, 1.0, 1.0), 100.0, 1.0, 0.0, 0.0, 0.0)},
	{"Red", materialStats(dvec3(1.0, 0.0, 0.0), 2.0, 1.0, 0.0, 0.0, 0.5)},
	{"Green", materialStats(dvec3(0.0, 1.0, 0.0), 100.0, 1.0, 0.0, 0.6, 0.0)},
	{"Bug", materialStats(dvec3(1.0, 0.0, 1.0), 100.0, 1.0, 0.0, 0.0, 1.0)},
	{"Copper", materialStats(dvec3(0.7038,0.27048,0.0828), 100.0, 1.0, 0.0, 0.5, 0.5)},
	{"Mirror", materialStats(dvec3(1.0, 1.0, 1.0), 50.0, 0.13511, 0.0, 1.0, 1.0)},
	{"MirrorB", materialStats(dvec3(0.7038,0.27048,0.0828), 50.0, 1.0, 0.0, 0.6, 0.6)},
};







//default bad constuctor
Material::Material() {
	this->color = dvec3(1.0, 0.0, 1.0);
	this->ns = 100.0;
	this->ni = 1.0;
	this->transparency = 0.0;
	this->metalness = 0.5;
	this->smoothness = 0.5;
}


Material::Material(
	dvec3 color,
	double ns,
	double ni,
	double transparency,
	double metalness,
	double smoothness,
	dvec3 emission,
	function<dvec3(dvec2)> colorFn,
	function<dvec3(dvec2)> normalFn,
	function<double(dvec2)> nsFn,
	function<double(dvec2)> niFn,
	function<double(dvec2)> transparencyFn,
	function<double(dvec2)> metalnessFn,
	function<double(dvec2)> smoothnessFn
) {
	this->color = color;
	this->ns = ns;
	this->ni = ni;
	this->transparency = transparency;
	this->metalness = metalness;
	this->smoothness = smoothness;
	this->emission = emission;

	this->colorFn = colorFn;
	this->normalFn = normalFn;
	this->niFn = niFn;
	this->nsFn = nsFn;
	this->transparencyFn = transparencyFn;
	this->metalnessFn = metalnessFn;
	this->smoothnessFn = smoothnessFn;
}
/*Material::Material(
	materialStats stats) {

	this->color = stats.color;
	this->ns = stats.ns;
	this->ni = stats.ni;
	this->transparency = stats.transparency;
	this->metalness = stats.metalness;
	this->smoothness = stats.smoothness;
}*/

Material::Material(
	string premade,
	dvec3 emmision,
	function<dvec3(dvec2)> colorFn,
	function<dvec3(dvec2)> normalFn,
	function<double(dvec2)> nsFn,
	function<double(dvec2)> niFn,
	function<double(dvec2)> transparencyFn,
	function<double(dvec2)> metalnessFn,
	function<double(dvec2)> smoothnessFn
) {
	materialStats stats = materials.at(premade);
	this->color = stats.color;
	this->ns = stats.ns;
	this->ni = stats.ni;
	this->transparency = stats.transparency;
	this->metalness = stats.metalness;
	this->smoothness = stats.smoothness;
	this->emission = emmision;
	this->colorFn = colorFn;
	this->normalFn = normalFn;
	this->niFn = niFn;
	this->nsFn = nsFn;
	this->transparencyFn = transparencyFn;
	this->metalnessFn = metalnessFn;
	this->smoothnessFn = smoothnessFn;
}



dvec3 Material::getColor(dvec2 uv) const {
	//printf("gettingColor\n");
	if (this->colorFn != nullptr) {
		//printf("colorfn exists\n");
		//printf("that fn's value: %s at uv %s\n", glm::to_string(this->colorFn(uv)).c_str(), glm::to_string(uv).c_str());
		return this->colorFn(uv);
	}
	//printf("just returning default color\n");
	return this->color;
}
dvec3 Material::getNormal(dvec2 uv, dvec3 defaultNormal) const {
	if (this->normalFn != nullptr) {
		return this->normalFn(uv);
	}
	return defaultNormal;
}
double Material::getNS(dvec2 uv) const {
	if (this->nsFn != nullptr) {
		return this->nsFn(uv);
	}
	return this->ns;

}
double Material::getNI(dvec2 uv) const {
	if (this->niFn != nullptr) {
		return this->niFn(uv);
	}
	return this->ni;
}
double Material::getTransparency(dvec2 uv) const {
	if (this->transparencyFn != nullptr) {
		return this->transparencyFn(uv);
	}
	return this->transparency;
}
double Material::getMetalness(dvec2 uv) const {
	if (this->metalnessFn != nullptr) {
		return this->metalnessFn(uv);
	}
	return this->metalness;
}
double Material::getSmoothness(dvec2 uv) const {
	if (this->smoothnessFn != nullptr) {
		return this->smoothnessFn(uv);
	}
	return this->smoothness;
}

dvec3 Material::getEmission() const {
	return this->emission;
}
Material* Material::setColor(dvec3 newColor) {
	this->color = newColor;
	return this;
}