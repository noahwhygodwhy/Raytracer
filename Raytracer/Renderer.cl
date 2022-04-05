#include "rand.cl"
#include "sharedStructs.cl"
#include "sphere.cl"
#include "shape.cl"
#include "helpers.cl"
#include "cookTorance.cl"



#define BIAS 1e-4f

HitResult rayHitListOfShapes(const Ray ray, __global const Sphere* spheres, uint numberOfSpheres){
    int x = (int)1<INFINITY;

    HitResult minRayResult;
	minRayResult.hit = false;;
	minRayResult.depth = INFINITY;

    for (uint i = 0; i < numberOfSpheres; i++) {
        
        Sphere theSphere = spheres[i];
        rayAABBResult resA = rayAABB(theSphere.shape.boundingBox, ray);

        if (resA.hit) {
            HitResult resB = rayHitSphere(theSphere, ray, i);
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


//TODO: this might not be an ok way of sampling hemisphere vectors
float3 randomHemisphericalVector(float3 normal, ulong* state) {



    //return (float3)(0.0, 1.0, 0.0);
    //return (float3)(0.0f, 1.0f, 0.0f);
//     float3 randomVector = normalize((float3)(rand(idx+randomCounter++), rand(idx+randomCounter++), rand(idx+randomCounter++)));
//     //if(idx == 500+(1000*550))printf("new random hemisphere vector: %f, %f, %f\n", randomVector.x, randomVector.y, randomVector.z);
//    // if(idx == 500+(1000*550))printf("normal vector: %f, %f, %f\n", normal.x, normal.y, normal.z);
//     float dotVN = dot(randomVector, normal);
//     if(fabs(dotVN)==1.0){
//         return (float3)(0.0, 1.0, 0.0);

//     } else {
//         float3 rotVec = cross(randomVector, normal);
//         return normalize(rotateVector(randomVector, acos(dotVN), rotVec));
        
//     }


    //return (float3)(0.0f, 1.0f, 0.0f);

    //printf("called random hemispherical vector\n");
	float3 outVec = (float3)(0.0f, 0.0f, 0.0f);

	normal = normalize(normal);
	do {
		outVec = normalize((float3)(randDubThree(state), randDubThree(state), randDubThree(state)));
	} while (dot(normal, outVec) < 0.0);
    //printf("random hemisphere vector from normal %f, %f, %f: %f, %f, %f\n", normal.x, normal.y, normal.z, outVec.x, outVec.y, outVec.z);
	return outVec;
}






HitResult shootRay(Ray ray, __global const Sphere* spheres, uint numberOfSpheres) {

	return rayHitListOfShapes(ray, spheres, numberOfSpheres);

    //kd tree stuff goes here if you do that
}

__kernel void render(
    __global const OtherData* otherData,
    __global const Sphere* spheres,
    __global const Material* materials,
    __global float4* frameBuffer,
    __global ulong* randomBuffer) //an array of maxJumps hitResults
    {

    float randomCounter = otherData->currentTime;

    int frameX = get_global_size(0);
    int frameY = get_global_size(1);

    int pixelX = get_global_id(0);
    int pixelY = get_global_id(1);


    int sampleN = 0;//get_global_id(2);


    int pixelIdx = pixelX+(frameX*pixelY);

    ulong state = randomBuffer[pixelIdx];

    // printf("state: %lu\n", state);

    float normalizedX = (((float)(pixelX)/(float)(frameX))-0.5f);//*2.0f;
    float normalizedY = (((float)(pixelY)/(float)(frameY))-0.5f);//*2.0f;


    float4 coordOnScreen = (normalizedX * otherData->camRight) + (normalizedY * otherData->camUp) + otherData->eye + (otherData->camForward * otherData->focal);
    
    float2 clipSpacePixelSize = (float2)(1.0f / max((float)(frameX - 1.0f), 1.0f), 1.0f / max((float)(frameY - 1.0f), 1.0f));


    float multiSampleX = 0.0f;//TODO:
    float multiSampleY = 0.0f;//TODO:


    float offsetX = (multiSampleX * ((clipSpacePixelSize.x / (float)(sampleN + 1))))+(clipSpacePixelSize.x/2.0f);//TODO: this is wrong
    float offsetY = (multiSampleY * ((clipSpacePixelSize.y / (float)(sampleN + 1))))+(clipSpacePixelSize.y/2.0f);

    Ray ray;
    Ray newRay;
// 
    ray.direction = (float3)normalize((coordOnScreen.xyz + (float3)(offsetX, offsetY, 0.0f) ) - otherData->eye.xyz);

    ray.origin = otherData->eye.xyz;
    newRay = ray;

//   printf("pixelidx: %i, randomseed: %lu, ^^^ %lu\n", pixelIdx, otherData->randomSeed, ((ulong)pixelIdx) ^ otherData->randomSeed); 
//     return;


//     printf("state: %u\n", otherData->randomSeed);

    
//     return;  
    
    //printf("state: %u\n", otherData->randomSeed);
    // return;
    // float3 totalRadiance = (float3)(0.0f, 0.0f, 0.0f);

    // float3 layerMultiplier = (float3)(1.0f, 1.0f, 1.0f);
    // float3 prevLayerRadiance = (float3)(0.0f, 0.0f, 0.0f);
    //float3 kS = (float3)(1.0f, 1.0f, 1.0f);
    //float3 cT = (float3)(1.0f, 1.0f, 1.0f);
   // float3 kD = (float3)(1.0f, 1.0f, 1.0f);

    //float3 upperColor = (float3)(1.0f, 1.0f, 1.0f);

    float ior = 1.0;

    //float3 thisRadiance =(float3)(0.0f, 0.0f, 0.0f);

    //HitResult prevHit;
    HitResult newHit;
    

    //http://raytracey.blogspot.com/2016/11/opencl-path-tracing-tutorial-2-path.html

    float3 monteAccum = (float3)(0.0f);
    //printf("number of samples: %u\n", otherData->numberOfSamples);

    //return;
    for(uint monte = 0; monte < otherData->numberOfSamples; monte++){
        ray.direction = (float3)normalize((coordOnScreen.xyz + (float3)(offsetX, offsetY, 0.0f) ) - otherData->eye.xyz);
        ray.origin = otherData->eye.xyz;
        newRay = ray;

        float3 accumulated = (float3)(0.0f);
        float3 masked = (float3)(1.0f);
        randomCounter+=1.5123f;
        for(uint layer = 0; layer < otherData->maxDepth; layer++) {

            ray = newRay;
            newHit = shootRay(ray, spheres, otherData->numberOfSpheres);

            if(!newHit.hit){
                if(layer==0u){
                    accumulated = masked*otherData->clearColor.xyz;
                }
                break;

            }

            float transparencyDecider = rand(&state);
            float reflectanceDecider = rand(&state);
            //float specularDecider = rand(pixelIdx+randomCounter++);

            float hitAngle = acos(dot(newHit.normal, -ray.direction));

            bool entering = hitAngle < (M_PI_F / 2.0);

            Material mat = materials[spheres[newHit.shapeIdx].shape.matIdx];


            //printf("spheres[newHit.shapeIdx].shape.matIdx %u, layer: %u, deciders: %f<%f, %f<%f, emmision: %f, %f, %f\n", spheres[newHit.shapeIdx].shape.matIdx, layer, transparencyDecider, mat.trans,  reflectanceDecider, mat.smooth, mat.emission.r, mat.emission.g, mat.emission.b);//, specularDecider);
            
            /*if(!(mat.emission.r == 0.0 && mat.emission.g == 0.0 && mat.emission.b == 0.0)){

                break;
            }*/
            
            if (transparencyDecider < mat.trans) {
                newRay.direction = normalize(getRefractionRay(normalize(newHit.normal), normalize(ray.direction), mat.ni, entering));
                newRay.origin = newHit.position + (newHit.normal * (entering ? -1.0f : 1.0f) * BIAS);
            }
            else {
                newRay.origin = newHit.position + (newHit.normal * BIAS);

                if (reflectanceDecider < mat.smooth && false) {
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

                // float3 F0a = (float3)(fabs((1.0f - mat.ni) / (1.0f + mat.ni)));
                // F0a = F0a * F0a;
                // float3 matC = mat.color.xyz;

                // float3 F0 = mix(F0a, matC, mat.metal);
                // //dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, minRayResult.shape->mat, currentIOR, F0, kS);
                // //TODO: need to keep the current ray
                // float3 kS;
                // float3 cT = CookTorance(ray, newRay, newHit, mat, ior, F0, &kS);

        //blinn phong (ish)
                //float3 R = rotateVector(newRay.direction, M_PI_F, newHit.normal);//TODO: not the right rotate
                float3 V = -ray.direction;
                float3 N = newHit.normal;
                float3 L = newRay.direction;

                float3 H = (L + V) / length(L + V);

                float diff = max(dot(L, N), 0.0f);
                float spec = max(dot(N, H), 0.0f);

                //float3 kD = ((1.0f - kS) * (1.0f - mat.metal))*diff;*/

                //if(pixelX == 825 && pixelY == 300)printf("V: %f, %f, %f\n", V.x, V.y, V.z);
                //if(pixelX == 825 && pixelY == 300)printf("L: %f, %f, %f\n", L.x, L.y, L.z);
                //if(pixelX == 825 && pixelY == 300)printf("N: %f, %f, %f\n", N.x, N.y, N.z);
                //if(pixelX == 825 && pixelY == 300)printf("diff: %f, spec: %f\n", diff, spec);

                accumulated += mat.emission.xyz*masked;
                masked *= mat.color.xyz;
                masked*= diff+spec;
                
                //accumulated = (float3)(fabs(newHit.normal.x), fabs(newHit.normal.y), fabs(newHit.normal.z));
                //break;
                // if(pixelIdx == 500+(1000*550))printf("KD, KS, CT: %f, %f, %f\n", kD.x, kS.y, cT.z);
                //if(pixelX == 500 && pixelY == 600)printf("%i, %i, KD, KS, CT: %f, %f, %f\n", pixelX, pixelY, kD.x, kS.y, cT.z);
                
                //totalRadiance = totalRadiance + (layerMultiplier * mat.color.xyz);
                //prevLayerRadiance = mat.color.xyz;
                //layerMultiplier *= (kD);

                
            }
            
            
        }
        //if(pixelX == 825 && pixelY == 300)printf("accumlated is %f, %f, %f\n", accumulated.x, accumulated.y, accumulated.z);
        monteAccum += accumulated;
    }
    //uint layer;

    //totalRadiance = (float3)(0.5, 0.0, 0.0);
    //float4 result = {1.0, 0.0, 0.0, 1.0};
    frameBuffer[pixelIdx] = (float4)((monteAccum)/(float)(otherData->numberOfSamples), 1.0f);
    randomBuffer[pixelIdx] = state;
    //if(pixelX == 825 || pixelY == 300) frameBuffer[pixelIdx] = (float4)(1.0, 0.0, 1.0, 1.0);
    //frameBuffer[pixelIdx] = (float4)(1.0, 0.0, 0.0, 1.0);
}


//slide 25 https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.2pp.pdf