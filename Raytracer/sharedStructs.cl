#ifndef SHARED_STRUCTS_H
#define SHARED_STRUCTS_H


#ifdef CPP

    typedef struct alignas(128) rayAABBResult {
        float4 enter;
        float4 exit;
        bool hit;
        
    }  rayAABBResult;



    typedef struct alignas(128) OtherData {
        float4 clearColor;
        float4 eye;
        float4 camRight;
        float4 camUp;
        float4 camForward;
        float focal;
        uint maxDepth;
        uint numberOfSpheres;
    }  OtherData;


    typedef struct alignas(128) Material {
        float4 color;
        float4 emission;
        float ns;
        float ni;
        float trans;
        float metal;
        float smooth;
    }  Material;



    typedef struct alignas(128) Ray {
        float4 origin;
        float4 direction;
    } Ray;


    typedef struct alignas(128) HitResult {
        float4 position;
        float4 normal;
        float4 uv;
        float depth;
        uint shapeIdx;
        bool hit;
    }  HitResult;

    typedef struct alignas(128) AABB {
        float4 min;
        float4 max; 
    }  AABB;


    typedef struct alignas(128) Shape {
        AABB boundingBox;
        uint matIdx;
        //uint shapeIdx;
        
    }  Shape;

    typedef struct alignas(128) Sphere {
        float4 origin;
        Shape shape;
        float radius;

        
    } Sphere;




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
        float4 enter;
        float4 exit;
        bool hit;
        
    } rayAABBResult;



    typedef struct OtherData {
        float4 clearColor;
        float4 eye;
        float4 camRight;
        float4 camUp;
        float4 camForward;
        float focal;
        uint maxDepth;
        uint numberOfSpheres;
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
        float4 origin;
        float4 direction;
    } Ray;


    typedef struct HitResult {
        float4 position;
        float4 normal;
        float4 uv;
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
#endif





#endif
