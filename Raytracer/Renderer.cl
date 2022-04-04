#include "rand.cl"
#include "sharedStructs.cl"
#include "sphere.cl"
#include "shape.cl"
#include "helpers.cl"
#include "cookTorance.cl"

//#define OFFSETOF(TYPE, ELEMENT) ((size_t)&(((TYPE *)0)->ELEMENT))

bool prd = false;


float bias = 1e-4;

HitResult rayHitListOfShapes(const Ray ray, const Sphere* spheres, uint numberOfSpheres){
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
float3 randomHemisphericalVector(float3 normal, int idx, int randomCounter) {

    //return (float3)(0.0, 1.0, 0.0);
    //return (float3)(0.0f, 1.0f, 0.0f);
    float3 randomVector = normalize((float3)(rand(idx+randomCounter+0), rand(idx+randomCounter+1), rand(idx+randomCounter+2)));
    //if(idx == 500+(1000*550))printf("new random hemisphere vector: %f, %f, %f\n", randomVector.x, randomVector.y, randomVector.z);
   // if(idx == 500+(1000*550))printf("normal vector: %f, %f, %f\n", normal.x, normal.y, normal.z);
    float dotVN = dot(randomVector, normal);
    if(fabs(dotVN)==1.0){
        return (float3)(0.0, 1.0, 0.0);

    } else {
        float3 rotVec = cross(randomVector, normal);
        return normalize(rotateVector(randomVector, acos(dotVN), rotVec));
        
    }


    //return (float3)(0.0f, 1.0f, 0.0f);

    //printf("called random hemispherical vector\n");
	// float3 outVec = (float3)(0.0f, 0.0f, 0.0f);

	// normal = normalize(normal);
	// do {
	// 	outVec = normalize((float3)(randDubThree(idx+randomCounter+0), randDubThree(idx+randomCounter+1), randDubThree(idx+randomCounter+2)));
	// } while (dot(normal, outVec) < 0.0);
    // //printf("random hemisphere vector from normal %f, %f, %f: %f, %f, %f\n", normal.x, normal.y, normal.z, outVec.x, outVec.y, outVec.z);
	// return outVec;
}






HitResult shootRay(Ray ray, const Sphere* spheres, uint numberOfSpheres) {

	return rayHitListOfShapes(ray, spheres, numberOfSpheres);

    //kd tree stuff goes here if you do that
}

__kernel void render(
    __global const OtherData* otherData,
    __global const Sphere* spheres,
    __global const Material* materials,
    __global float4* frameBuffer) //an array of maxJumps hitResults
    {

    float randomCounter = otherData->currentTime;

    int frameX = get_global_size(0);
    int frameY = get_global_size(1);

    int pixelX = get_global_id(0);
    int pixelY = get_global_id(1);


    int sampleN = 0;//get_global_id(2);

    prd = (pixelX == 500) && (pixelY == 500);

    int pixelIdx = pixelX+(frameX*pixelY);


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


    

    float3 totalRadiance = (float3)(0.0f, 0.0f, 0.0f);

    float3 layerMultiplier = (float3)(1.0f, 1.0f, 1.0f);
    float3 prevLayerRadiance = (float3)(0.0f, 0.0f, 0.0f);
    //float3 kS = (float3)(1.0f, 1.0f, 1.0f);
    //float3 cT = (float3)(1.0f, 1.0f, 1.0f);
   // float3 kD = (float3)(1.0f, 1.0f, 1.0f);

    //float3 upperColor = (float3)(1.0f, 1.0f, 1.0f);

    float ior = 1.0;

    //float3 thisRadiance =(float3)(0.0f, 0.0f, 0.0f);

    //HitResult prevHit;
    HitResult newHit;
    

    float3 accumulated = (float3)(0.0f);
    float3 masked = (float3)(1.0f);

    
    for(uint layer = 0; layer < otherData->maxDepth; layer++) {
        
        ray = newRay;
        newHit = shootRay(ray, spheres, otherData->numberOfSpheres);

        if(!newHit.hit){
            if(layer==1u){
                accumulated = masked*otherData->clearColor.xyz;
            }
            break;

        }

        float transparencyDecider = rand(pixelIdx+randomCounter+1);
        float reflectanceDecider = rand(pixelIdx+randomCounter+1);
        //float specularDecider = rand(pixelIdx+randomCounter++);

        float hitAngle = acos(dot(newHit.normal, -ray.direction));

        bool entering = hitAngle < (M_PI_F / 2.0);

        Material mat = materials[spheres[newHit.shapeIdx].shape.matIdx];


        //printf("spheres[newHit.shapeIdx].shape.matIdx %u, layer: %u, deciders: %f<%f, %f<%f, emmision: %f, %f, %f\n", spheres[newHit.shapeIdx].shape.matIdx, layer, transparencyDecider, mat.trans,  reflectanceDecider, mat.smooth, mat.emission.r, mat.emission.g, mat.emission.b);//, specularDecider);
        ;
        /*if(!(mat.emission.r == 0.0 && mat.emission.g == 0.0 && mat.emission.b == 0.0)){

            break;
        }*/
        
        if (transparencyDecider < mat.trans) {
            newRay.direction = normalize(getRefractionRay(normalize(newHit.normal), normalize(ray.direction), mat.ni, entering));
            newRay.origin = newHit.position + (newHit.normal * (entering ? -1.0f : 1.0f) * bias);
        }
        else {

            if (reflectanceDecider < mat.smooth) {
                newRay.direction = -rotateVector(ray.direction, M_PI_F, newHit.normal);//this is not the correct "rotate"
                newRay.origin = newHit.position + (newHit.normal * bias);
            }
            else {
                
                newRay.direction = normalize(randomHemisphericalVector(newHit.normal, pixelIdx, randomCounter));
                randomCounter+=3;
                newRay.origin = newHit.position + (newHit.normal * bias);
            }

 
//calculate lighting stuff

            float3 F0a = (float3)(fabs((1.0f - mat.ni) / (1.0f + mat.ni)));
            F0a = F0a * F0a;
            float3 matC = mat.color.xyz;

            float3 F0 = mix(F0a, matC, mat.metal);
            //dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, minRayResult.shape->mat, currentIOR, F0, kS);
            //TODO: need to keep the current ray
            float3 kS;
            float3 cT = CookTorance(ray, newRay, newHit, mat, ior, F0, &kS);

    //blinn phong (ish)
            float3 R = rotateVector(newRay.direction, M_PI_F, newHit.normal);//TODO: not the right rotate
            float3 V = ray.origin;
            float3 N = newHit.normal;
            float3 L = newRay.direction;

            float3 H = (L + V) / length(L + V);

            float diff = dot(L, N);
            float spec = dot(N, H);

            float3 kD = ((1.0f - kS) * (1.0f - mat.metal))*diff;


            accumulated += mat.emission.xyz*masked;
            masked *= mat.color.xyz *(diff*kD);
            //accumulated = (float3)(diff*kD);
            break;
            // if(pixelIdx == 500+(1000*550))printf("KD, KS, CT: %f, %f, %f\n", kD.x, kS.y, cT.z);
            //if(pixelX == 500 && pixelY == 600)printf("%i, %i, KD, KS, CT: %f, %f, %f\n", pixelX, pixelY, kD.x, kS.y, cT.z);
            
            //totalRadiance = totalRadiance + (layerMultiplier * mat.color.xyz);
            //prevLayerRadiance = mat.color.xyz;
            //layerMultiplier *= (kD);

            
        }
        
    }

    //uint layer;

    //totalRadiance = (float3)(0.5, 0.0, 0.0);
    //float4 result = {1.0, 0.0, 0.0, 1.0};
    frameBuffer[pixelIdx] = (float4)(accumulated, 1.0)/otherData->numberOfSamples;
    //if(pixelX == 500 && pixelY == 600) frameBuffer[pixelIdx] = (float4)(1.0, 0.0, 0.0, 1.0);
}


//slide 25 https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.2pp.pdf