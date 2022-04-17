#include "rand.cl"
#include "sharedStructs.cl"
#include "triangle.cl"
#include "sphere.cl"
#include "shape.cl"
#include "helpers.cl"
#include "cookTorance.cl"
#include "material.cl"




#define VOLUME_SAMPLE_DISTANCE 0.05f

#define BIAS 1e-3f



HitResult rayHitListOfShapes(const Ray ray, 
uint numberOfShapes,
__global const UShape* shapes, 
__global const Vertex* vertices){
    int x = (int)1<INFINITY;

    HitResult minRayResult;
	minRayResult.hit = false;;
	minRayResult.depth = INFINITY;

    for(uint i = 0; i < numberOfShapes; i++){
        rayAABBResult resA = rayAABB(shapes[i].boundingBox, ray);
        if (resA.hit) {
            HitResult resB;
            
            switch(shapes[i].type){
            case 0:
                resB = rayHitSphere(shapes[i], ray, i);
                break;
            case 1:
                resB = rayHitTriangle(shapes[i], ray, vertices, i);
                break;
            default:
                return minRayResult;
                break;
            }
            
            if(resB.hit) {
                if (resB.depth < minRayResult.depth) {
                    minRayResult = resB;
                }
            }
        

        }
    }
    return minRayResult;
}    


float3 getRefractionRay(float3 hitNormal, float3 incidentVector, float objectIOR, bool entering) {


	float closeness = dot(hitNormal.xyz, incidentVector.xyz);
	float prevIOR = 1.0;
	float newIOR = objectIOR;


	if (!entering) {
		hitNormal = -hitNormal;
        float s = prevIOR;
        prevIOR = newIOR;
        newIOR = s;
	}

	float cosA1 = dot(incidentVector, hitNormal);

	float sinA1 = sqrt(1.0 - (cosA1 * cosA1));

	float IORRatio = prevIOR / newIOR;

	float sinA2 = sinA1 * IORRatio;

	float3 trueReflectDir = incidentVector;
	if (sinA2 <= -1.0 || sinA2 >= 1.0) {//TODO this isn't handled correctly
		return trueReflectDir;
	}
	float maxCloseness = -INFINITY;
	float k1 = NAN;
	float k2 = NAN;
	float2 temp = solveQuadratic(1.0, 2.0*cosA1, 1.0-(1.0/(IORRatio * IORRatio)));
    k1 = temp.x;
    k2 = temp.y;
    
	if (!isnan(k1)) {
		float3 reflectDir1 = normalize(incidentVector + (k1 * hitNormal));
		float closeness1 = dot(incidentVector, reflectDir1);
		if (closeness1 > maxCloseness && closeness1>=0.0) {
			maxCloseness = closeness1;
			trueReflectDir = reflectDir1;
		}
	}
	if (!isnan(k2)) {
		float3 reflectDir2 = normalize(incidentVector + (k2 * hitNormal));
		float closeness2 = dot(incidentVector, reflectDir2);
		if (closeness2 > maxCloseness && closeness2>=0.0) {
			maxCloseness = closeness2;
			trueReflectDir = reflectDir2;
		}
	}

	if (maxCloseness <= 0.0) {
		printf("error calculating refraction angle\n");
		return incidentVector;
		//exit(-1);
	}

	float cosA2 = sqrt(1.0 - (sinA2 * sinA2));

	if (cosA1 < 0.0) {
		cosA2 *= -1.0;
	}

	return normalize(trueReflectDir);
}


float3 randomSphericalVector(ulong* state){
    return normalize((float3)(randDubThree(state), randDubThree(state), randDubThree(state)));
}
//TODO: this might not be an ok way of sampling hemisphere vectors
float3 randomHemisphericalVector(float3 normal, ulong* state) {


	float3 outVec = (float3)(0.0f, 0.0f, 0.0f);

	normal = normalize(normal);
	do {
		outVec = normalize((float3)(randDubThree(state), randDubThree(state), randDubThree(state)));
	} while (dot(normal, outVec) < 0.0);
    //printf("random hemisphere vector from normal %f, %f, %f: %f, %f, %f\n", normal.x, normal.y, normal.z, outVec.x, outVec.y, outVec.z);
	return outVec;
}






HitResult shootRay(Ray ray, uint numberOfShapes,
__global const UShape* shapes, __global const Vertex* vertices) {
	return rayHitListOfShapes(ray, numberOfShapes, shapes, vertices);
    //kd tree stuff goes here if you do that
}

__kernel void render(
    __global const OtherData* otherData,
    __global const UShape* shapes,
    __global const Vertex* vertices,
    __global const Material* materials,
    __global ulong* randomBuffer,
    read_write image2d_t outBuffer,
    __global volatile ToneMapStruct* toneMapInfo,
    const uint frameCounter
    //__global float2* minMax
    )  {




    float currentFrame = ((float)frameCounter)/(float)otherData->fps;

    int frameX = get_global_size(0);
    int frameY = get_global_size(1);

    float frameRatio = (float)frameX/(float)frameY;
    float3 eye;
    eye = (float3)(sin(currentFrame) * 40, 12, cos(currentFrame) * 40);

    //eye = (float3)(0.0f, 7.0f, 40.0f);
    float3 lookat = (float3)(0.0, 0.0, 0.0);

    float3 camForward = normalize(lookat - eye);
    float3 camUp = normalize((float3)(0.0, 1, 0.0));
    float3 camRight = cross(camForward, camUp);
    camUp = cross(camRight, camForward);

    float viewPortHeight = 2.0f;
    float viewPortWidth = viewPortHeight * frameRatio;

    float fov = 90.0f;
    float focal = (viewPortHeight / 2.0) / tan(radians(fov / 2.0));


    int pixelX = get_global_id(0);
    int pixelY = get_global_id(1);

    int sampleN = 0;//get_global_id(2);


    int pixelIdx = pixelX+(frameX*pixelY);
    


    ulong state = randomBuffer[pixelIdx];


    float normalizedX = (((float)(pixelX)/(float)(frameX))-0.5f);//*2.0f;
    float normalizedY = (((float)(pixelY)/(float)(frameY))-0.5f);//*2.0f;


    float4 coordOnScreen = (float4)((normalizedX * camRight) + (normalizedY * camUp) + eye + (camForward * focal), 0.0);
    
    float2 clipSpacePixelSize = (float2)(1.0f / max((float)(frameX - 1.0f), 1.0f), 1.0f / max((float)(frameY - 1.0f), 1.0f));


    float multiSampleX = 0.0f;//TODO:
    float multiSampleY = 0.0f;//TODO:


    float offsetX = (multiSampleX * ((clipSpacePixelSize.x / (float)(sampleN + 1))))+(clipSpacePixelSize.x/2.0f);//TODO: this is wrong
    float offsetY = (multiSampleY * ((clipSpacePixelSize.y / (float)(sampleN + 1))))+(clipSpacePixelSize.y/2.0f);

    Ray ray;
    Ray newRay;
// 
    // ray.direction = (float3)normalize((coordOnScreen.xyz + (float3)(offsetX, offsetY, 0.0f) ) - eye.xyz);

    // ray.origin = eye.xyz;
    // newRay = ray;


    float ior = 1.0;

    HitResult newHit;
    

    //http://raytracey.blogspot.com/2016/11/opencl-path-tracing-tutorial-2-path.html

    float3 monteAccum = (float3)(0.0f);
    for(uint monte = 0; monte < otherData->numberOfSamples; monte++){
        ray.direction = (float3)normalize((coordOnScreen.xyz + (float3)(offsetX, offsetY, 0.0f) ) - eye.xyz);
        ray.origin = eye.xyz;
        newRay = ray;

        float3 accumulated = (float3)(0.0f);
        float3 masked = (float3)(1.0f);
        bool insideFog = false;

        for(uint layer = 0; layer < otherData->maxDepth; layer++) {

            if(false)printf("\nlayer %u\n", layer);
            ray = newRay;
            
            newHit = shootRay(ray, otherData->numberOfShapes, shapes, vertices);
            //if(pixelX == 500 && pixelY == 500)printf("hit: %i\n", (int)newHit.hit);
            //hit needsmat idx not shape idx
            if(!newHit.hit){
                if(false)printf("Didn't hit nuffin\n");
                if(layer == 0u){
                    accumulated = masked*otherData->clearColor.xyz;
                }
                break;

            }
            
            if(false)printf("hit something\n");
            // frameBuffer[pixelIdx] = (float4)(1.0, 0.0, 1.0, 1.0);
            // return;
            //printf("there was a hit\n");
            
            float transparencyDecider = rand(&state);
            float reflectanceDecider = rand(&state);


            float hitAngle = acos(dot(newHit.normal, -ray.direction));

            bool entering = hitAngle < (M_PI_F / 2.0);


            Material mat = materials[newHit.matIdx];


            float3 matColor = getColor(mat.color.xyz, newHit.uv.xy, mat.proceduralColor, currentFrame);








            if(mat.volumetricColor > 0u && entering){
                
                if(false)printf("hit fog case 1\n");
                insideFog = true;
                //if(pixelX == 500 && pixelY == 500) printf("entering\n");
                newRay.origin = newHit.position + (-newHit.normal * BIAS);
                newRay.direction = ray.direction;
            } 
            else if(insideFog){
                
                if(false)printf("already inside fog\n");
                newRay.origin = newHit.position+ (-newHit.normal * BIAS);
                newRay.direction = ray.direction;
                
                if(false)printf("max stepper distance: %f\n", newHit.depth);
                for(float stepper = 0.0f; stepper < newHit.depth; stepper+=VOLUME_SAMPLE_DISTANCE){
                    if(false)printf("layer %u, stepper at %f\n", layer, stepper);

                    float scatterDecider = rand(&state);
                    
                    //if(scatterDecider < mat.color.w) {
                    if(scatterDecider < getFogDensity(mat.volumetricColor, ray.origin+(ray.direction*stepper) + (newHit.normal * BIAS), currentFrame, mat.color.w)) {
                        if(false)printf("scattering\n" );
                        newRay.origin = ray.origin+(ray.direction*stepper) + (newHit.normal * BIAS);
                        newRay.direction = randomSphericalVector(&state);
                        break;
                    }
                    
                    newRay.origin = newHit.position+ (newHit.normal * BIAS);
                    newRay.direction = ray.direction;
                    insideFog = false;
                }
            }



//TODO these three cases
/*outside fog and entering
inside fog but just bounding
inside fog and exiting*/

            // if(mat.volumetricColor > 0u || insideFog){
            //     if(entering && !insideFog){
            //         insideFog = true;
            //         //if(pixelX == 500 && pixelY == 500) printf("entering\n");
            //         newRay.origin = newHit.position + (-newHit.normal * BIAS);
            //         newRay.direction = ray.direction;
            //     } else {
                    
            //         newRay.origin = newHit.position + (newHit.normal * BIAS);
            //         newRay.direction = ray.direction;
            //         //if(pixelX == 500 && pixelY == 500) printf("exiting, going to depth %f\n", newHit.depth);

            //         for(float stepper = 0.0f; stepper < newHit.depth; stepper+=VOLUME_SAMPLE_DISTANCE){

            //             float scatterDecider = rand(&state);
            //             //if(pixelX == 500 && pixelY == 500)printf("stepper: %f, colorw: %f, scatterDecider: %f\n", stepper, mat.color.w, scatterDecider);
                        
            //             if(scatterDecider < getFogDensity(mat.volumetricColor, newHit.position, currentFrame, mat.color.w)) {
            //                 //if(pixelX == 500 && pixelY == 500)printf("scattering\n");
            //                 newRay.origin = ray.origin+(ray.direction*stepper) + (newHit.normal * BIAS);
            //                 newRay.direction = randomSphericalVector(&state);
            //                 break;
            //             }
            //         }
            //         insideFog = false;
            //     }

            // }


            else if (transparencyDecider < mat.trans) {
                
                if(false)printf("transparency\n");
                newRay.direction = normalize(getRefractionRay(normalize(newHit.normal), normalize(ray.direction), mat.ni, entering));
                newRay.origin = newHit.position + (newHit.normal * (entering ? -1.0f : 1.0f) * BIAS);
            }
            else {
                
                if(false)printf("not transparent\n");
                newRay.origin = newHit.position + (newHit.normal * BIAS);

                if (reflectanceDecider < mat.smooth) {
                    newRay.direction = -rotateVector(ray.direction, M_PI_F, newHit.normal);//this is not the correct "rotate"
                }
                else {
                    
                    
                    //TODO: this was essentially stolen from 
                    //http://raytracey.blogspot.com/2016/11/opencl-path-tracing-tutorial-2-path.html
                    //need to go through and see how the one I wrote was wrong after I confirm this works.
                    float3 normal = newHit.normal;
                    float3 normal_facing = dot(normal, ray.direction) < 0.0f ? normal : normal * (-1.0f);

                    /* compute two random numbers to pick a random point on the hemisphere above the hitpoint*/
                    
                    float rand1 = 2.0f * M_PI_F * rand(&state);
                    float rand2 = rand(&state);
                    float rand2s = sqrt(rand2);

                    /* create a local orthogonal coordinate frame centered at the hitpoint */
                    float3 w = normal_facing;
                    float3 axis = fabs(w.x) > 0.1f ? (float3)(0.0f, 1.0f, 0.0f) : (float3)(1.0f, 0.0f, 0.0f);
                    float3 u = normalize(cross(axis, w));
                    float3 v = cross(w, u);
                    
                    float3 newdir = normalize(u * cos(rand1)*rand2s + v*sin(rand1)*rand2s + w*sqrt(1.0f - rand2));
                    
                    newRay.direction = newdir;//normalize(randomHemisphericalVector(newHit.normal, &state));
                }


    
    //calculate lighting stuff

                float3 F0a = (float3)(fabs((1.0f - mat.ni) / (1.0f + mat.ni)));
                F0a = F0a * F0a;

                float3 F0 = mix(F0a, matColor, mat.metal);
                // //dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, minRayResult.shape->mat, currentIOR, F0, kS);
                // //TODO: need to keep the current ray
                float3 kS;
                float3 cT = max(CookTorance(ray, newRay, newHit, mat, ior, F0, &kS), 0.0f);
        //blinn phong (ish)
                //float3 R = rotateVector(newRay.direction, M_PI_F, newHit.normal);//TODO: not the right rotate
                float3 V = -ray.direction;
                float3 N = newHit.normal;
                float3 L = newRay.direction;

                float3 H = normalize((L + V));

                float diff = max(dot(L, N), 0.0f);
                float spec = max(dot(N, H), 0.0f);

                float3 kD = ((1.0f - kS) * (1.0f - mat.metal));//*diff;
                
                accumulated += mat.emission.xyz*masked;
                masked *= matColor;
                masked*= (diff*kD)+(cT);  
            }
            
            

            // if(mat.volumetricColor > 0u && !entering){
                
            //     if(pixelX == 600 || pixelY == 800)printf("leaving the fog?\n");
            //     newRay.origin = newHit.position + (newHit.normal * BIAS);
            //     newRay.direction = ray.direction;
            //     insideFog = false;
            // }
            
        }

        monteAccum += accumulated;
    }
    //printf("on pixel %i\n", pixelIdx);
    float3 result = (monteAccum)/(float)(otherData->numberOfSamples);


    
    //if(pixelX == 600 || pixelY == 800) result = (float3)(1.0f, 0.0f, 1.0f);

    write_imagef(outBuffer, (int2) (pixelX, pixelY),  (float4)(result, 1.0f));

    struct onionRGB onionResults;
    onionResults.r.f = result.x;
    onionResults.g.f = result.y;
    onionResults.b.f = result.z;
    
    atomic_min(&(toneMapInfo->minLum.r.i), onionResults.r.i);
    atomic_min(&(toneMapInfo->minLum.g.i), onionResults.g.i);
    atomic_min(&(toneMapInfo->minLum.b.i), onionResults.b.i);

    atomic_max(&(toneMapInfo->maxLum.r.i), onionResults.r.i);
    atomic_max(&(toneMapInfo->maxLum.g.i), onionResults.g.i);
    atomic_max(&(toneMapInfo->maxLum.b.i), onionResults.b.i);
    
    struct onionRGB logOnionResults;

    float oneOverN = 1.0f/frameX*frameY;


    float absoluteLum = 0.27*result.x + 0.67*result.y * 0.06 * result.z;


    logOnionResults.r.f = log(0.0001+absoluteLum) * oneOverN;
    logOnionResults.g.f = log(0.0001+absoluteLum) * oneOverN;
    logOnionResults.b.f = log(0.0001+absoluteLum) * oneOverN;
    
    atomic_add(&(toneMapInfo->summedLogLum.r.i), logOnionResults.r.i);
    atomic_add(&(toneMapInfo->summedLogLum.g.i), logOnionResults.g.i);
    atomic_add(&(toneMapInfo->summedLogLum.b.i), logOnionResults.b.i);

    randomBuffer[pixelIdx] = state;
}






__kernel void toneMap(
    __global const OtherData* otherData,
    read_write image2d_t outBuffer,
    __global volatile ToneMapStruct* toneMapInfo
    )  {



    int frameX = get_global_size(0);
    int frameY = get_global_size(1);

    int pixelX = get_global_id(0);
    int pixelY = get_global_id(1);

    int sampleN = 0;//get_global_id(2);

    int pixelIdx = pixelX+(frameX*pixelY);




//ward tone reproduction method
    float3 maxes = (float3) (toneMapInfo->maxLum.r.f, toneMapInfo->maxLum.g.f, toneMapInfo->maxLum.b.f);
    float3 logavg = exp((float3) (toneMapInfo->summedLogLum.r.f, toneMapInfo->summedLogLum.g.f, toneMapInfo->summedLogLum.b.f));
    float3 pix = read_imagef(outBuffer, (int2) (pixelX, pixelY)).xyz;
    float3 top = 1.219f +pow(maxes/2.0f, 0.4f);
    float3 bot = 1.219f  +pow(logavg, 0.4f);
    float3 result = pix * (1.0f/maxes)*pow(top/bot, 2.5f);


//reinhard tone reproduction method


    // float3 pix = read_imagef(outBuffer, (int2) (pixelX, pixelY)).xyz;

    // float absoluteLum = 0.27*pix.x + 0.67*pix.y * 0.06 * pix.z;


    // float LwLocal = exp(log(0.0001 + absoluteLum));



    // float3 maxes = (float3) (toneMapInfo->maxLum.r.f, toneMapInfo->maxLum.g.f, toneMapInfo->maxLum.b.f);
    // float3 LwTotal = exp((float3)(toneMapInfo->summedLogLum.r.f, toneMapInfo->summedLogLum.g.f, toneMapInfo->summedLogLum.b.f));
    // float a = 0.18f;
    // float3 scaled = (a/LwTotal)*LwLocal;
    // float3 reflectance = scaled/(1.0f+scaled);
    // float3 result = reflectance*maxes;


    // if(pixelIdx == 0){

    //     printf("maxes: %f, %f, %f, maxes 2: %f, logavg: %f,%f,%f, logavg2: %f", maxes.x, maxes.y, maxes.z, maxes2, logavg.x, logavg.y, logavg.z, logavg2);

    // }

    //result = pix;

    result = pow(result, (float3)(1.0f/2.2f));//gamma correction for sRGB
    write_imagef(outBuffer, (int2) (pixelX, pixelY), (float4)(result, 1.0f));

}


//slide 25 https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.2pp.pdf