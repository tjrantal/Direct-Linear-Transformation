package timo.home.dlt;

import static org.jocl.CL.*;

import java.util.Arrays;

import org.jocl.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;

public class JOCL_views_to_3D {

    //views_to_global(const int camNo,const int points, __global const float* coordinates, __global const float* coefficients, __global float* reco)
    public static float[] testJOCL(int camNo, int points, float[] coordinates, float[] coefficients){
        String programSource = "";
        float[] reco = new float[points*3]; //return array
        try{
            programSource = readKernelFile("kernels/views_to_3d.cl");
            //System.out.println(programSource);
        }catch(Exception e){
            System.out.println("Spat the dummy");
            return null;
        }


        Pointer pcoordinates = Pointer.to(coordinates);
        Pointer pcoefficients = Pointer.to(coefficients);
        Pointer preco = Pointer.to(reco);

        // The platform, device type and device number
        // that will be used
        final int platformIndex = 0;
        final long deviceType = CL_DEVICE_TYPE_GPU;
        final int deviceIndex = 0;

        // Enable exceptions and subsequently omit error checks in this sample
        CL.setExceptionsEnabled(true);

        // Obtain the number of platforms
        int numPlatformsArray[] = new int[1];
        clGetPlatformIDs(0, null, numPlatformsArray);
        int numPlatforms = numPlatformsArray[0];

        // Obtain a platform ID
        cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
        clGetPlatformIDs(platforms.length, platforms, null);
        cl_platform_id platform = platforms[platformIndex];

        // Initialize the context properties
        cl_context_properties contextProperties = new cl_context_properties();
        contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);
        
        // Obtain the number of devices for the platform
        int numDevicesArray[] = new int[1];
        clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
        int numDevices = numDevicesArray[0];
        
        // Obtain a device ID 
        cl_device_id devices[] = new cl_device_id[numDevices];
        clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
        cl_device_id device = devices[deviceIndex];

        // Create a context for the selected device
        cl_context context = clCreateContext(
            contextProperties, 1, new cl_device_id[]{device}, 
            null, null, null);
        
        // Create a command-queue for the selected device
        cl_queue_properties properties = new cl_queue_properties();
        cl_command_queue commandQueue = clCreateCommandQueueWithProperties(
            context, device, properties, null);

        // Allocate the memory objects for the input- and output data
        cl_mem srcCoordinates = clCreateBuffer(context, 
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            Sizeof.cl_float * coordinates.length, pcoordinates, null);
        cl_mem srcCoefficients = clCreateBuffer(context, 
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            Sizeof.cl_float * coefficients.length, pcoefficients, null);
        cl_mem dstReco = clCreateBuffer(context, 
            CL_MEM_READ_WRITE, 
            Sizeof.cl_float * reco.length, null, null);
        
        // Create the program from the source code
        cl_program program = clCreateProgramWithSource(context,
            1, new String[]{ programSource }, null, null);
        
        // Build the program
        clBuildProgram(program, 0, null, null, null, null);
        
        // Create the kernel
        cl_kernel kernel = clCreateKernel(program, "views_to_global", null);
        
        // Set the arguments for the kernel. Use pointers to convey the required arguments
        int a = 0;
        clSetKernelArg(kernel, a++, Sizeof.cl_int, Pointer.to(new int[]{camNo}));
        clSetKernelArg(kernel, a++, Sizeof.cl_int, Pointer.to(new int[]{points}));
        clSetKernelArg(kernel, a++, Sizeof.cl_mem, Pointer.to(srcCoordinates));
        clSetKernelArg(kernel, a++, Sizeof.cl_mem, Pointer.to(srcCoefficients));
        clSetKernelArg(kernel, a++, Sizeof.cl_mem, Pointer.to(dstReco));
        
        // Set the work-item dimensions
        long global_work_size[] = new long[]{points};
        
        // Execute the kernel
        clEnqueueNDRangeKernel(commandQueue, kernel, global_work_size.length, null,
            global_work_size, null, 0, null, null);
        
        // Read the output data into reco
        clEnqueueReadBuffer(commandQueue, dstReco, CL_TRUE, 0,
            reco.length*Sizeof.cl_float, preco, 0, null, null);
        
        // Release kernel, program, and memory objects
        clReleaseMemObject(srcCoordinates);
        clReleaseMemObject(srcCoefficients);
        clReleaseMemObject(dstReco);
        clReleaseKernel(kernel);
        clReleaseProgram(program);
        clReleaseCommandQueue(commandQueue);
        clReleaseContext(context);
        
        /* 
            // Verify the result
            float[] check = mat_mul(n,srcArrayA,srcArrayB);
            boolean passed = true;
            final float epsilon = 1e-7f;
            for (int i=0; i<srcArrayA.length; i++)
            {
                float x = dstArray[i];
                float y = check[i];
                boolean epsilonEqual = Math.abs(x - y) <= epsilon * Math.abs(x);
                if (!epsilonEqual)
                {
                    passed = false;
                    break;
                }
            }
            System.out.println("Test "+(passed?"PASSED":"FAILED"));
            if (n <= srcArrayA.length)
            {
                System.out.println("GPU: "+Arrays.toString(dstArray));
                System.out.println("CPU: "+Arrays.toString(check));
            }
        */
        return reco;
    }

    private static float[] mat_mul(int N,float[] a, float[] b){
        float[] c = new float[a.length];
        int i, j, k;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                c[i*N+j] = 0.0f;
                for (k = 0; k < N; k++) {
                    // C(i, j) = sum(over k) A(i,k) * B(k,j)
                    c[i*N+j] += a[i*N+k] * b[k*N+j];
                }
            }
        }
        return c;
    }

    public static String readKernelFile(String path) throws IOException {
        ClassLoader classLoader = Thread.currentThread().getContextClassLoader();
        try (InputStream inputStream = classLoader.getResourceAsStream(path);
             InputStreamReader streamReader = new InputStreamReader(inputStream, StandardCharsets.UTF_8);
             BufferedReader reader = new BufferedReader(streamReader)) {

            StringBuilder stringBuilder = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                stringBuilder.append(line).append("\n");
            }
            return stringBuilder.toString();
        }
    }
}
