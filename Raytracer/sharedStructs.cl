#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H


enum shapeType {
	SPHERE = 0,
	TRIANGLE


};

#ifdef CPP

    

    typedef struct alignas(16) OtherData {
        
	    float4 clearColor;
        //float4 eye;
        //float4 camRight;
        //float4 camUp;
        //float4 camForward;
        //float focal;
        //float currentTime;
        uint fps;
        uint maxDepth;
        uint numberOfShapes;
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
        uint proceduralColor;
        uint volumetricColor;
    }  Material;

    typedef struct alignas(16) AABB {
        float4 min;
        float4 max; 
    }  AABB;

    typedef struct alignas(16) UShape{
        float4 values;
        AABB boundingBox;
        uint matIdx;
        uint type;
    } UShape;

    typedef struct alignas(16) Shape {
        AABB boundingBox;
        uint matIdx;
        uint type;
        
    }  Shape;

    // typedef struct alignas(16) Sphere {
    //     float4 origin;
    //     Shape shape;
    //     float radius;

        
    // } Sphere;

    typedef struct Vertex{
        float4 position;
        float4 normal;
        float4 uv;
    } Vertex;

    // typedef struct Triangle{
    //     Shape shape;
    //     uint vertA;
    //     uint vertB;
    //     uint vertC;
    // }Triangle;



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
        //float4 eye;
        //float4 camRight;
        //float4 camUp;
        //float4 camForward;
        //float focal;
        //float currentTime;
        uint fps;
        uint maxDepth;
        uint numberOfShapes;
        uint numberOfSamples;
    }  OtherData;


    typedef struct Material {
        float4 color;
        float4 emission;
        float ns;
        float ni;
        float trans;
        float metal;
        float smooth;
        uint proceduralColor;
        uint volumetricColor;
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
        uint matIdx;
        bool hit;
    } HitResult;

    typedef struct AABB {
        float4 min;
        float4 max; 
    } AABB;


    typedef struct UShape{
        float4 values;
        AABB boundingBox;
        uint matIdx;
        uint type;
    } UShape;
    typedef struct Shape {
        AABB boundingBox;
        uint matIdx;
        uint type;
        
    } Shape;

    // typedef struct Sphere {
    //     float4 origin;
    //     Shape shape;
    //     float radius;

        
    // } Sphere;

    typedef struct Vertex{
        float4 position;
        float4 normal;
        float4 uv;
    } Vertex;

    // typedef struct Triangle{
    //     Shape shape;
    //     uint vertA;
    //     uint vertB;
    //     uint vertC;
    // }Triangle;
#endif





#endif
