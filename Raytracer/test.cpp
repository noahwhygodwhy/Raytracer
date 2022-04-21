//
//#include <cstdint>
//#include <glm/glm.hpp>
//using namespace glm;
//typedef vec4 float4;
//
//
//union atomFloat {
//    float f;
//    int i;
//};
//
//typedef struct alignas(16) OtherData {
//
//    float4 clearColor;
//    uint fps;
//    uint maxDepth;
//    uint numberOfShapes;
//    uint numberOfSamples;
//}  OtherData;
//
//
//
//typedef struct alignas(16) Material {
//    float4 color;
//    float4 emission;
//    float ns;
//    float ni;
//    float trans;
//    float metal;
//    float smooth;
//    uint proceduralColor;
//    uint volumetricColor;
//}  Material;
//
//typedef struct alignas(16) AABB {
//    float4 min;
//    float4 max;
//}  AABB;
//
//typedef struct alignas(16) GKDNode
//{
//    AABB area;
//    float value;
//    int greater;
//    int lesser;
//    int axis;
//} GKDNode;
//
//typedef struct alignas(16) UShape {
//    float4 values;
//    AABB boundingBox;
//    uint matIdx;
//    uint type;
//    int next;
//} UShape;
//
//
//typedef struct alignas(16) Vertex {
//    float4 position;
//    float4 normal;
//    float4 uv;
//} Vertex;
//
//typedef struct alignas(16) onionRGB {
//    union atomFloat r;
//    union atomFloat g;
//    union atomFloat b;
//
//} onionRGB;
//
//typedef struct alignas(16) ToneMapStruct {
//    struct onionRGB minLum;
//    struct onionRGB maxLum;
//    struct onionRGB summedLogLum;
//}ToneMapStruct;
