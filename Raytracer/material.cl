#ifndef MATERIAL_H
#define MATERIAL_H


float3 getColor(float3 colorIn, float2 uv, float swc) {
    
    switch((int)floor(swc)) {
        case 0:
            return colorIn;
            break;
        case 1:
            return (float3)(uv.x, uv.y, 1.0f);
            break;
        case 2:
            //float2 squares(10.0, 10.0);
            // int2 flatUV = (int2)floor(uv * 10);
            if (fmod(((floor(uv.x*20.0f))+ (floor(uv.y*30.0f))),2) == 0) {
                return (float3)(1.0, 0.0, 0.0);

            }
            else {
                return (float3)(1.0, 1.0, 0.0);
            }
            return (float3)(uv.x, uv.y, 1.0f);
            break;
        default:
            return (float3)(1.0f, 0.0f, 1.0f);
    }
}

#endif
