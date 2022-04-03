

__kernel void square(__global const int* a, __global int* b) {
    int gid = get_global_id(0);
    b[gid] = (a[gid]*a[gid]);
}


//slide 25 https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.2pp.pdf