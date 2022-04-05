#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H



#ifdef CPP

    //typedef float4 cl_float4;
    //typedef float cl_float;
    // typedef struct alignas(16) rayAABBResult {
    //     float4 enter;
    //     float4 exit;
    //     bool hit;
        
    // }  rayAABBResult;



    typedef struct alignas(16) OtherData {
        float4 clearColor;
        float4 eye;
        float4 camRight;
        float4 camUp;
        float4 camForward;
        cl_long randomSeed;
        float focal;
        float currentTime;
        uint maxDepth;
        uint numberOfSpheres;
        uint numberOfTriangles;
        uint numberOfSamples;
    }  OtherData;


    typedef struct alignas(16) Material {
        float4 color;
        float4 emission;
        float ns;
        float ni;
        float trans;
        float metal;
        float smooth;
    }  Material;



    // typedef struct alignas(16) Ray {
    //     float4 origin;
    //     float4 direction;
    // } Ray;


    // typedef struct alignas(16) HitResult {
    //     float4 position;
    //     float4 normal;
    //     float4 uv;
    //     float depth;
    //     uint shapeIdx;
    //     bool hit;
    // }  HitResult;

    typedef struct alignas(16) AABB {
        float4 min;
        float4 max; 
    }  AABB;


    typedef struct alignas(16) Shape {
        AABB boundingBox;
        uint matIdx;
        //uint shapeIdx;
        
    }  Shape;

    typedef struct alignas(16) Sphere {
        float4 origin;
        Shape shape;
        float radius;

        
    } Sphere;

    typedef struct Vertex{
        float4 position;
        float4 normal;
        float4 uv;
    } Vertex;

    typedef struct Triangle{
        uint vertA;
        uint vertB;
        uint vertC;
    }Triangle;



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
#else
    
    

    typedef struct rayAABBResult {
        float3 enter;
        float3 exit;
        bool hit;
        
    } rayAABBResult;



    typedef struct OtherData {
        
	    float4 clearColor;
        float4 eye;
        float4 camRight;
        float4 camUp;
        float4 camForward;
        ulong randomSeed;
        float focal;
        float currentTime;
        uint maxDepth;
        uint numberOfSpheres;
        uint numberOfTriangles;
        uint numberOfSamples;
    } OtherData;


    typedef struct Material {
        float4 color;
        float4 emission;
        float ns;
        float ni;
        float trans;
        float metal;
        float smooth;
    } Material;



    typedef struct Ray {
        float3 origin;
        float3 direction;
    } Ray;


    typedef struct HitResult {
        float3 position;
        float3 normal;
        float2 uv;
        float depth;
        uint shapeIdx;
        bool hit;
    } HitResult;

    typedef struct AABB {
        float4 min;
        float4 max; 
    } AABB;


    typedef struct Shape {
        AABB boundingBox;
        uint matIdx;
        //uint shapeIdx;
        
    } Shape;

    typedef struct Sphere {
        float4 origin;
        Shape shape;
        float radius;

        
    } Sphere;

    typedef struct Vertex{
        float4 position;
        float4 normal;
        float4 uv;
    } Vertex;

    typedef struct Triangle{
        uint vertA;
        uint vertB;
        uint vertC;
    }Triangle;
#endif





#endif
