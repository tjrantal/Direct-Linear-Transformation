__kernel void mat_mul(const int N, __global float *A,
__global float *B, __global float *C){
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k;
    float tmp = 0.0f;
    for (k = 0; k < N; k++) {
        // C(i, j) = sum(over k) A(i,k) * B(k,j)
        tmp += A[i*N+k] * B[k*N+j];
    }
    C[i*N+j] = tmp;
} 