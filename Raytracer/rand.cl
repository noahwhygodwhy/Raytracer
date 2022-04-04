




float randwithfloat2(float2 co){

    float2 t = {12.9898, 78.233};
    float a = dot(co, t);
    float b = sin(a);
    float c;
    return modf(b * 43758.5453, &c);
}


float rand(int xid){
    float a = ((((xid+5)<<6) * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1));
    float b = ((((xid+6)<<2) * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1));
    float2 c = (float2)(a, b);
    return fabs(randwithfloat2(c));
}

float randDubThree(int xid) {
    return (rand(xid)*2.0f)-1.0f;

}