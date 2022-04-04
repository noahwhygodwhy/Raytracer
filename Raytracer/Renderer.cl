#include "rand.cl"
#include "sharedStructs.cl"
#include "sphere.cl"
#include "shape.cl"
#include "helpers.cl"
#include "cookTorance.cl"



bool prd = false;


//float bias = 1e-4;


HitResult rayHitListOfShapes(const Ray ray, const Sphere* spheres, uint numberOfSpheres){
    HitResult minRayResult;
	minRayResult.hit = false;;
	minRayResult.depth = INFINITY;

    for (uint i = 0; i < numberOfSpheres; i++) {
        HitResult rayResult;
        
        Sphere theSphere = spheres[i];

        rayAABBResult resA = rayAABB(theSphere.shape.boundingBox, ray);

        if (resA.hit) {
            HitResult resB = rayHitSphere(theSphere, ray, i);
            if(resB.hit) {
                if (resB.depth < minRayResult.depth) {
                    minRayResult = rayResult;
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
	float3 outVec = (float3)(0.0f, 0.0f, 0.0f);
	normal = normalize(normal);
	do {
		outVec = normalize((float3)(randDubThree(idx+randomCounter+0), randDubThree(idx+randomCounter+1), randDubThree(idx+randomCounter+2)));
	} while (dot(normal, outVec) < 0.0);
	return outVec;
}



const char* printv(float4 a) {
    //char* x;
    printf("%f, %f, %f, %f\n", a.x, a.y, a.z, a.w);
    //return x;

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

    int randomCounter = 0;

    int frameX = get_global_size(0);
    int frameY = get_global_size(1);

    int pixelX = get_global_id(0);
    int pixelY = get_global_id(1);
    int sampleN = 0;//get_global_id(2);

    int pixelIdx = pixelX+(frameX*pixelY);


    float normalizedX = ((float)(pixelX)/(float)(frameX))-0.5f;
    float normalizedY = ((float)(pixelY)/(float)(frameY))-0.5f;

    float4 coordOnScreen = (normalizedX * otherData->camRight) + (normalizedY * otherData->camUp) + otherData->eye + (otherData->camForward * otherData->focal);
    float2 clipSpacePixelSize = (float2)(1.0f / (float)(frameX - 1.0f), 1.0f / (float)(frameY - 1.0f));

    float offsetX = pixelX * (clipSpacePixelSize.x / (sampleN + 1));
    float offsetY = pixelY * (clipSpacePixelSize.y / (sampleN + 1));

    //float4 result = (float4)(rand(pixelIdx),rand(pixelIdx),rand(pixelIdx),rand(pixelIdx));

    Ray ray;
    Ray newRay;

    ray.direction = (float4)(normalize((coordOnScreen.xyz + (float3)(offsetX, offsetY, 0.0f)) - otherData->eye.xyz), 0.0);
    ray.origin = otherData->eye;

    
    

    float3 totalRadiance = (float3)(0.0f, 0.0f, 0.0f);
    float3 layerMultiplier = (float3)(1.0f, 1.0f, 1.0f);
    float3 kS = (float3)(1.0f, 1.0f, 1.0f);
    float3 cT = (float3)(1.0f, 1.0f, 1.0f);
    float3 kD = (float3)(1.0f, 1.0f, 1.0f);
    float3 upperColor = (float3)(1.0f, 1.0f, 1.0f);

    float ior = 1.0;

    float3 thisRadiance =(float3)(0.0f, 0.0f, 0.0f);

    HitResult prevHit;
    HitResult newHit;
    
    newHit = shootRay(ray, spheres, otherData->numberOfSpheres);

    //printf("%f, %f, %f\n", otherData->eye.x, otherData->eye.y, otherData->eye.z);
    Sphere x = spheres[0];

    prd = pixelX == 500 && pixelY == 500;
    if(prd) {
        printf("1\n");
        printf("2\n");

        printf("origin:");
        printv(x.origin);
        printv(x.shape.boundingBox.min);
        printv(x.shape.boundingBox.max);
        printf("matidx: %u\n", x.shape.matIdx);
        printf("radius: %u\n", x.radius);

    }
    // if("prd:")
    // if(prd)printf("origin:");
    // if(prd)printv(x.origin);
    // if(prd)printv(x.shape.boundingBox.min);
    // if(prd)printv(x.shape.boundingBox.max);
    // if(prd)printf("matidx: %u\n", x.shape.matIdx);
    // if(prd)printf("radius: %u\n", x.radius);
    
    if(newHit.hit) {

        //Material mat = materials[spheres[newHit.shapeIdx].shape.matIdx];
        frameBuffer[pixelIdx] = (float4)(1.0f, 0.0f, 0.0f, 1.0f);
    } else {
        frameBuffer[pixelIdx] = (float4)(0.0f, 0.0f, 1.0f, 1.0f);
    }
    /*if(otherData->numberOfSpheres > 0) {

        frameBuffer[pixelIdx] = (float4)(1.0, 0.0, 0.0, 1.0f);

    } else {

        frameBuffer[pixelIdx] = (float4)(0.0, 0.0f, 1.0f, 1.0f);

    }*/
    return;


    //frameBuffer[pixelIdx] = (float4)(1.0f, 1.0f, 0.0f, 1.0f);
    //return;
    for(uint layer = 0; layer < otherData->maxDepth; layer++) {


        //ray = newRay;
        newHit = shootRay(ray, spheres, otherData->numberOfSpheres);

        if(newHit.hit) {

            Material mat = materials[spheres[newHit.shapeIdx].shape.matIdx];
            frameBuffer[pixelIdx] = (float4)(mat.color.xyz, 1.0f);
        } else {
            frameBuffer[pixelIdx] = (float4)(otherData->clearColor.x, 1.0f, 0.0f, 1.0f);
        }
        return;
        //break;
        
        /*

            float transparencyDecider = rand(pixelIdx+randomCounter++);
            float reflectanceDecider = rand(pixelIdx+randomCounter++);
            float specularDecider = rand(pixelIdx+randomCounter++);
            float hitAngle = acos(dot(hits[layer].normal, -ray.direction));

            bool entering = hitAngle < (M_PI_F / 2.0);

            Material mat = materials[spheres[hits[layer].shapeIdx].shape.matIdx];

        }*/










/*

        if(hits[layer].hit) {


            float transparencyDecider = rand(pixelIdx+randomCounter++);
            float reflectanceDecider = rand(pixelIdx+randomCounter++);
            float specularDecider = rand(pixelIdx+randomCounter++);

            float hitAngle = acos(dot(hits[layer].normal, -ray.direction));

            bool entering = hitAngle < (M_PI_F / 2.0);

            Material mat = materials[spheres[hits[layer].shapeIdx].shape.matIdx];

            if(mat.emission.r == 0.0 && mat.emission.g == 0.0 && mat.emission.b == 0.0){
                if (transparencyDecider < mat.trans) {

                    newRay.direction = normalize(getRefractionRay(normalize(hits[layer].normal), normalize(ray.direction), mat.ni, entering));
                    newRay.origin = hits[layer].position + (hits[layer].normal * (entering ? -1.0f : 1.0f) * bias);

                    //thisRadiance = pathTrace(newRay, mat.ni, layer + 1); //TODO:

                }
                else {

                    if (reflectanceDecider < mat.smooth) {
                        newRay.direction = -rotateVector(ray.direction, M_PI_F, hits[layer].normal);//this is not the correct "rotate"
                        newRay.origin = hits[layer].position + (hits[layer].normal * bias);
                    }
                    else {
                        newRay.direction = normalize(randomHemisphericalVector(hits[layer].normal, pixelIdx, randomCounter));
                        randomCounter+=3;
                        newRay.origin = hits[layer].position + (hits[layer].normal * bias);
                    }

                    //downstreamRadiance = pathTrace(newRay, mat.ni, layer + 1); //TODO:



                    float3 kS = (float3)(0.0f, 0.0f, 0.0f);

                    float3 F0a = (float3)(fabs((1.0f - mat.ni) / (1.0f + mat.ni)));
                    F0a = F0a * F0a;
                    float3 matC = mat.color;

                    float3 F0 = mix(F0a, matC, mat.metal);
                    //dvec3 ctSpecular = CookTorance(ray, reflectionRay, minRayResult, downstreamRadiance, minRayResult.shape->mat, currentIOR, F0, kS);
                    //TODO: need to keep the current ray
                    cT = CookTorance(ray, newRay, hits[layer], mat, ior, F0, &kS);





            //blinn phong (ish)
                    float3 R = rotateVector(newRay.direction, M_PI_F, hits[layer].normal);//TODO: not the right rotate
                    float3 V = ray.origin;
                    float3 N = hits[layer].normal;
                    float3 L = newRay.direction;

                    float3 H = (L + V) / length(L + V);

                    float diff = dot(L, N);
                    float spec = dot(N, H);

                    float3 kD = (1.0f - kS) * (1.0f - mat.metal);;

                    upperColor = mat.color;
                }
            } else {
                
            }

        } else {
            totalRadiance += otherData->clearColor*layerMultiplier;
        }*/
    }

    //uint layer;

    
    //float4 result = {1.0, 0.0, 0.0, 1.0};
    //frameBuffer[gIdx] = result;
}


//slide 25 https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.2pp.pdf