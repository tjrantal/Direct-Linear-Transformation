#define maxMatrixRows 2*36 //camNo*2 needed at minimum, allow max 36 cams;

/*
    Solve group of equations using pseudo inverse (implementation from ChatGPT)
    matrixRows = rows in matrix A (and vector y)
    A = matrix A (matrixRows x width)
    y = solution to Ax = y
    x = the return value 
*/
void pseudo_inverse(int matrixRows, float* A, float* y, float* x) {
    int i, j, k;

    // Step 1: Compute A^T (transpose of A)
    float AT[3][maxMatrixRows];
    for (i = 0; i < matrixRows; ++i) {
        for (j = 0; j < 3; ++j) {
            AT[j][i] = A[i * 3 + j];
        }
    }

    // Step 2: Compute AT * A (resulting in 3x3 matrix)
    float ATA[3][3] = {0};
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < matrixRows; ++k) {
                ATA[i][j] += AT[i][k] * A[k * 3 + j];
            }
        }
    }

    // Step 3: Compute the inverse of ATA
    // Load ATA into local variables for convenience
    float a00 = ATA[0][0], a01 = ATA[0][1], a02 = ATA[0][2];
    float a10 = ATA[1][0], a11 = ATA[1][1], a12 = ATA[1][2];
    float a20 = ATA[2][0], a21 = ATA[2][1], a22 = ATA[2][2];

    // Compute the determinant of ATA
    float det = a00 * (a11 * a22 - a21 * a12) 
              - a01 * (a10 * a22 - a20 * a12) 
              + a02 * (a10 * a21 - a20 * a11);

    // If the determinant is zero, the matrix is not invertible.
    if (det == 0.0f) {
        // Handle non-invertible case (could set to identity matrix or handle error)
        // For now, we'll just set it to the identity matrix
        ATA[0][0] = 1.0f; ATA[0][1] = 0.0f; ATA[0][2] = 0.0f;
        ATA[1][0] = 0.0f; ATA[1][1] = 1.0f; ATA[1][2] = 0.0f;
        ATA[2][0] = 0.0f; ATA[2][1] = 0.0f; ATA[2][2] = 1.0f;
    } else {
        // Compute the adjugate of ATA (transposed cofactor matrix)
        float adj[3][3];
        adj[0][0] =  (a11 * a22 - a21 * a12);
        adj[0][1] = -(a01 * a22 - a21 * a02);
        adj[0][2] =  (a01 * a12 - a11 * a02);
        adj[1][0] = -(a10 * a22 - a20 * a12);
        adj[1][1] =  (a00 * a22 - a20 * a02);
        adj[1][2] = -(a00 * a12 - a10 * a02);
        adj[2][0] =  (a10 * a21 - a20 * a11);
        adj[2][1] = -(a00 * a21 - a20 * a01);
        adj[2][2] =  (a00 * a11 - a10 * a01);

        // Divide each element of the adjugate by the determinant to get the inverse
        float inv_det = 1.0f / det;
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                ATA[i][j] = adj[i][j] * inv_det;
            }
        }
    }

    // Step 4: Compute AT * y (resulting in 3x1 vector)
    float ATy[3] = {0};
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < matrixRows; ++j) {
            ATy[i] += AT[i][j] * y[j];
        }
    }

    // Step 5: Compute x = ATA_inv * ATy
    for (i = 0; i < 3; ++i) {
        x[i] = 0;
        for (j = 0; j < 3; ++j) {
            x[i] += ATA[i][j] * ATy[j];
        }
    }
}

/*
    Kernel to run through all points and camera pairs. Gid goes through each point.
    Can feed in 0 for missing data as long as at least 2 camera views have valid data
    camNo = number of camera views
    points = number of points to consider
    coordinates = pointer to camera coordinate values. interleaved cam0 x, cam0 y, cam1 x, cam1 y etc
*/
__kernel void views_to_global(const int camNo,const int points, __global const float* coordinates, __global const float* coefficients, __global float* reco){
    int i,pointStride,camStride,recoStride,coeffStride;
    int gid = get_global_id(0); //Get index of execution
    float L1[3*maxMatrixRows]  = {0};    //A into the pseudo_inverse
    float L2[maxMatrixRows]  = {0};  //y into the pseudo_inverse
    float x[3] = {0};    //Return the reconstructed coordinates
    pointStride = gid*camNo*2;  // 2 coordinates per view x views x current point
    //Prep matrices for solving Ax = y
    for (i =0;i<camNo;++i){
        camStride = i*2;
        coeffStride = i*11;
        L2[2*i] = coefficients[coeffStride+3]- coordinates[pointStride+camStride+0];
        L2[2*i+1] = coefficients[coeffStride+7]- coordinates[pointStride+camStride+1];
    }
    
    
    for (int i = 0;i<camNo;++i){
        camStride = i*2;
        coeffStride = i*11;
        L1[2*3*i+0]	=coefficients[coeffStride+8]*coordinates[pointStride+camStride+0]-coefficients[coeffStride+0];
        L1[2*3*i+1]	=coefficients[coeffStride+9]*coordinates[pointStride+camStride+0]-coefficients[coeffStride+1];
        L1[2*3*i+2]	=coefficients[coeffStride+10]*coordinates[pointStride+camStride+0]-coefficients[coeffStride+2];
        L1[2*3*i+3]	=coefficients[coeffStride+8]*coordinates[pointStride+camStride+1]-coefficients[coeffStride+4];
        L1[2*3*i+4]	=coefficients[coeffStride+9]*coordinates[pointStride+camStride+1]-coefficients[coeffStride+5];
        L1[2*3*i+5]	=coefficients[coeffStride+10]*coordinates[pointStride+camStride+1]-coefficients[coeffStride+6];
    }
    
    pseudo_inverse(camNo*2,L1,L2,x);   //Solve the group of equations
    recoStride = gid*3;
    for (i =  0;i<3;++i){
        reco[recoStride+i]=x[i]; //set the return values
    }
}
