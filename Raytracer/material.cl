#ifndef MATERIAL_H
#define MATERIAL_H


float3 getColor(float3 colorIn, float2 uv, float swc, float currentFrame) {
    
    //float fTime = fract(currentFrame);
    float biggerPiece;

    switch((int)floor(swc)) {
        case 0:
            return colorIn;
            break;
        case 1:
            return (float3)(uv.x, uv.y, fract(currentFrame, &biggerPiece));
            break;
        case 2:
            //float2 squares(10.0, 10.0);
            // int2 flatUV = (int2)floor(uv * 10);
            if (fmod(((floor(uv.x*10.0f))+ (floor(uv.y*10.0f))),2) == 0) {
                return (float3)(1, 0, 0);
            }
            else {
                return (float3)(1.0, 1.0,0);
            }
            break;
        case 3:
            //float2 squares(10.0, 10.0);
            // int2 flatUV = (int2)floor(uv * 10);
            if (fmod(((floor(uv.x*10.0f))+ (floor(uv.y*10.0f))),2) == 0) {
                return (float3)(0.2, 0.2, 0.2);
            }
            else {
                return (float3)(uv.x, uv.y, fract(currentFrame, &biggerPiece));
            }
            break;
        default:
            return (float3)(1.0f, 0.0f, 1.0f);
    }
}

float getFogDensity(float swc, float3 hitPos, float currentFrame, float originalDensity){

    float3 tempSphereOrigin = (float3)(7.5f, 6.0f, 7.5f);


    float sphereRad = 4.0f;

    switch((int)floor(swc)){
       
        case 3:{
        
            //return originalDensity;
            //return originalDensity;
            float wave1 = (sin(sin(hitPos.x) + (hitPos.y*cos(hitPos.z)) + currentFrame)+1.0f)/2.0f;
            float wave2 = 1.0f;
            float wave3 = 1.0f;
            //float wave1 = sin((hitPos.y*sin(((hitPos.z*hitPos.x)+currentFrame)*3.0f/sphereRad+currentFrame)+currentFrame)*5.0f/sphereRad)+0.1f;
            //float wave2 = sin((hitPos.z*sin(((hitPos.y)+currentFrame)*5.0f/sphereRad+currentFrame)+currentFrame)*3.0f/sphereRad)+0.2f;
            //float wave3 = cos((hitPos.z*cos(((hitPos.y)+currentFrame)*5.0f/sphereRad+currentFrame)+currentFrame)*3.0f/sphereRad)+0.2f;
            return wave1*originalDensity;
        }
        case 2:{
            
            if(fmod(((floor(hitPos.x/3.0f))+ (floor(hitPos.y/3.0f))+(floor(hitPos.z/3.0f))),2) == 0){
                return 0;
            }
            return originalDensity;
        }
        default:{
            return originalDensity;
        }
    }
}

#endif
