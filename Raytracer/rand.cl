
float randwithfloat2(float2 co){

    float2 t = {12.9898, 78.233};
    float a = dot(co, t);
    float b = sin(a);
    float c;
    return modf(b * 43758.5453, &c);
}


uint MWC64X(ulong *state)
{
    uint c=(*state)>>32, x=(*state)&0xFFFFFFFF;
    *state = x*((ulong)4294883355U) + c;
    return x^c;
}


float rand(ulong *state){
    
    float x = MWC64X(state)/((float)(UINT_MAX));
    //printf("new rand: %f\n", x);
    return x;
    // float a = ((((xid+5)<<6) * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1));
    // float b = ((((xid+6)<<2) * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1));
    // float2 c = (float2)(a, b);
    // return fabs(randwithfloat2(c));
}

float randDubThree(ulong* state) {
    return (rand(state)*2.0f)-1.0f;

}