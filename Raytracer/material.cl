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
    switch((int)floor(swc)){
        case 2:

            //return originalDensity*fabs(hitPos.y/10.0f);

            if(fmod(((floor(hitPos.x/3.0f))+ (floor(hitPos.y/3.0f))+(floor(hitPos.z/3.0f))),2) == 0){
                return 0;
            }
            return originalDensity;
        break;
        case 0://continue
        case 1://continue
        default:
            return originalDensity;
        break;
    }
}

#endif
