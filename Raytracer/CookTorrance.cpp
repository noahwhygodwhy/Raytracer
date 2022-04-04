//
//#include "CookTorrance.hpp"
//
//
//using namespace std;
//using namespace glm;
//
//
//extern bool prd;
//
//double Xchi(double x) {
//	int i = int(x > 0);
//	return (1.0 * double(i)) + (-1.0 * (double(1 - i)));
//
//}
//
//double D(dvec3 m, dvec3 n, dvec3 h, double a) {
//	double numerator = a * a * Xchi(glm::dot(h, n));
//	constexpr double pie = glm::pi<double>();
//
//	double mndot2 = glm::pow(glm::dot(m, n), 2);
//
//	double otherPart = (1.0 - mndot2) / mndot2;
//	double denominator = pie * glm::pow(glm::pow(glm::dot(m, n), 2)*((a*a)+otherPart), 2);
//
//	return numerator / denominator;
//}
//
//
//double fDialectric(double prevIOR, double newIOR, dvec3 i, dvec3 o) {
//	double f0 = (prevIOR - newIOR) / (prevIOR + newIOR);
//	f0 *= f0;
//
//	if (prd)printf("f0: %f, right: %f\n", f0, ((1.0 - f0) * glm::pow(1.0 - glm::dot(i, o), 5)));
//
//	return f0 + ((1.0-f0)*glm::pow(1.0 - glm::dot(i, o), 5));
//
//}
//dvec3 fresnelSchlick(double cosT, dvec3 f0)
//{
//	return f0 + (1.0 - f0) * pow(1.0 - cosT, 5);
//}
//
//double k = 4;
//
//double fConductor(double matIOR, dvec3 i, dvec3 o) {
//	double numerator = glm::pow(matIOR - 1.0, 2) + (4.0 * matIOR * glm::pow(1.0 - glm::dot(i, o), 5)) + glm::pow(k, 2);
//	double denom = glm::pow(matIOR + 1.0, 2) + glm::pow(k, 2);
//
//	if (prd)printf("num: %f, denom: %f\n", numerator, denom);
//	return numerator / denom;
//
//}
//
//
//double Gp(dvec3 v, dvec3 n, dvec3 h, double a) {
//
//
//	double vdh = dot(v, h);
//	double chi = Xchi(vdh / dot(v, n));
//	vdh *= vdh;
//	double right = (1.0 - vdh) / vdh;
//	return (chi * 2.0) / (1.0 + sqrt(1.0 + (a * a * right)));
//}
//double G(dvec3 wi, dvec3 wo, dvec3 n, dvec3 h, double a) {
//	return Gp(wi, n, h, a) * Gp(wo, n, h, a);
//}
//
//
//dvec3 CookTorance(
//	const Ray& incomingRay,
//	const Ray& outgoingRay,
//	const HitResult& minRayResult,
//	dvec3 downstreamRadiance,
//	Material* mat,
//	double prevIOR,
//	const dvec3& F0,
//	dvec3& kS) {
//
//	if (prd)printf("ct\n");
//
//	dvec2 uv = minRayResult.uv;
//	
//	double roughness = 1.0 - mat->getSmoothness(uv);
//
//
//	dvec3 vi = incomingRay.inverseDirection;
//	dvec3 vo = outgoingRay.direction;
//	//dvec3 vh = glm::normalize(vi + vo);
//	dvec3 vn = minRayResult.normal;
//
//	/*if (glm::dot(vn, vo) < 0) { //if going into transparent material
//		vo = -vo;
//		vo = vo - ((2.0 * vo * vn) / glm::pow(glm::length(vo), 2) * vo);
//	} else if(glm::dot(vn, vi) < 0) { //if coming out of transparent material
//		vi = -vi;
//		vi = vi - ((2.0 * vi * vn) / glm::pow(glm::length(vi), 2) * vi);
//	}*/
//	dvec3 vh = glm::normalize(vi + vo);
//
//	if (prd)printf("vi: %s\n", glm::to_string(vi).c_str());
//	if (prd)printf("vo: %s\n", glm::to_string(vo).c_str());
//	if (prd)printf("vh: %s\n", glm::to_string(vh).c_str());
//	if (prd)printf("vn: %s\n", glm::to_string(vn).c_str());
//
//	double ndi = glm::dot(vn, vi);
//
//	double cosT = glm::dot(vo, vn);
//	double sinT = glm::sqrt(1 - cosT * cosT);
//
//	// Calculate fresnel
//	dvec3 fresnel = fresnelSchlick(glm::dot(vh, vi), F0);
//
//	if(prd)printf("fresnel: %s\n", glm::to_string(fresnel).c_str());
//
//	// Geometry term
//	double geometry = G(vi, vo, vn, vh, roughness);// GGX_PartialGeometryTerm(viewVector, normal, halfVector, roughness)* GGX_PartialGeometryTerm(sampleVector, normal, halfVector, roughness);
//	if (prd)printf("geometry: %f\n", geometry);
//												   
//												   // Calculate the Cook-Torrance denominator
//	double denominator = 4.0 * (ndi * glm::dot(vh, vn) + 0.05);
//	if (prd)printf("denominator: %f\n", denominator);
//	kS = fresnel;
//	// Accumulate the radiance
//	return downstreamRadiance * geometry * fresnel * sinT / denominator;
//}