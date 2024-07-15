/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

    Copyright (C) 2012-2024 Timo Rantalainen
*/
/*
	3D DLT (direct linear transformation).
	Depends on Apache math https://commons.apache.org/proper/commons-math/
	Written by Timo Rantalainen tjrantal@gmail.com
*/

package timo.home.dlt;

//Apache Math
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;


public class DLT3D_apache_math3{
	
	public static double[] getDltCoefficients(double[][] calibrationObjectGlobalCoordinates,double[][] L){
		//Check for missing coordinates
		int validCoords = 0;
		for (int i = 0; i<L.length;++i){
			if (L[i] != null){
				++validCoords;
			}
		}
		//Indices for valid coordinates
		int[] validIndices = new int[validCoords];
		validCoords = 0;
		for (int i = 0; i<L.length;++i){
			if (L[i] != null){
				validIndices[validCoords++] = i;
			}	
		}
		
	
	
		double[][] B = new double[2*validIndices.length][11];//Matrix for solving DLT-parameters
		double[] C = new double[2*validIndices.length]; //Digitized calibrationObject coordinates
		int cnt = 0;
		for (int j =0;j<validIndices.length;j++){
			for (int i =0;i<2;i++){
				C[cnt++] = L[validIndices[j]][i];
			}
		}
	
		for (int i=0;i<validIndices.length;i++){
			B[2*i][0]	=calibrationObjectGlobalCoordinates[validIndices[i]][0];
			B[2*i][1]	=calibrationObjectGlobalCoordinates[validIndices[i]][1];
			B[2*i][2]	=calibrationObjectGlobalCoordinates[validIndices[i]][2];
			B[2*i][3]	=1;
			B[2*i][4]	=0;
			B[2*i][5]	=0;
			B[2*i][6]	=0;
			B[2*i][7]	=0;
			B[2*i][8]	=-calibrationObjectGlobalCoordinates[validIndices[i]][0]*L[validIndices[i]][0];
			B[2*i][9]	=-calibrationObjectGlobalCoordinates[validIndices[i]][1]*L[validIndices[i]][0];
			B[2*i][10]	=-calibrationObjectGlobalCoordinates[validIndices[i]][2]*L[validIndices[i]][0];
			B[2*i+1][0]	=0;
			B[2*i+1][1]	=0;
			B[2*i+1][2]	=0;
			B[2*i+1][3]	=0;
			B[2*i+1][4]	=calibrationObjectGlobalCoordinates[validIndices[i]][0];
			B[2*i+1][5]	=calibrationObjectGlobalCoordinates[validIndices[i]][1];
			B[2*i+1][6]	=calibrationObjectGlobalCoordinates[validIndices[i]][2];
			B[2*i+1][7]	=1;
			B[2*i+1][8]	=-calibrationObjectGlobalCoordinates[validIndices[i]][0]*L[validIndices[i]][1];
			B[2*i+1][9]	=-calibrationObjectGlobalCoordinates[validIndices[i]][1]*L[validIndices[i]][1];
			B[2*i+1][10]=-calibrationObjectGlobalCoordinates[validIndices[i]][2]*L[validIndices[i]][1];
		}
		return (solveVectorUsingSVD(MatrixUtils.createRealMatrix(B), MatrixUtils.createRealVector(C))).toArray();
	}

	public static double[] scaleCoordinates(double[][] coefficients, double[][] coordinates){
		double[][]	L1 = new double[2*coefficients.length][3];
		double[]		L2 = new double[2*coordinates.length];
		
		
		for (int i =0;i<coordinates.length;++i){
			L2[2*i] = coefficients[i][3]- coordinates[i][0];
			L2[2*i+1] = coefficients[i][7]- coordinates[i][1];
		}
		
		for (int i = 0;i<coefficients.length;++i){
			L1[2*i][0]	=coefficients[i][8]*coordinates[i][0]-coefficients[i][0];
			L1[2*i][1]	=coefficients[i][9]*coordinates[i][0]-coefficients[i][1];
			L1[2*i][2]	=coefficients[i][10]*coordinates[i][0]-coefficients[i][2];
			L1[2*i+1][0]	=coefficients[i][8]*coordinates[i][1]-coefficients[i][4];
			L1[2*i+1][1]	=coefficients[i][9]*coordinates[i][1]-coefficients[i][5];
			L1[2*i+1][2]	=coefficients[i][10]*coordinates[i][1]-coefficients[i][6];
		}
		return (solveVectorUsingSVD(MatrixUtils.createRealMatrix(L1), MatrixUtils.createRealVector(L2))).toArray();
	}

	//Project known 3D coordinates into the camera view
	public static double[] project3D(double[] coefficients, double[] globalCoordinates){
		double u = (coefficients[0]*globalCoordinates[0]+coefficients[1]*globalCoordinates[1]+coefficients[2]*globalCoordinates[2]+coefficients[3])/(coefficients[8]*globalCoordinates[0]+coefficients[9]*globalCoordinates[1]+coefficients[10]*globalCoordinates[2]+1);
		double v = (coefficients[4]*globalCoordinates[0]+coefficients[5]*globalCoordinates[1]+coefficients[6]*globalCoordinates[2]+coefficients[7])/(coefficients[8]*globalCoordinates[0]+coefficients[9]*globalCoordinates[1]+coefficients[10]*globalCoordinates[2]+1);
		return new double[]{u,v};
	}

	//Get camera position from DLT coefficients
	public static double[] getCameraPosition(double[] coefficients){
		double[][] A = new double[][] {{coefficients[0], coefficients[1], coefficients[2]},{coefficients[4], coefficients[5], coefficients[6]} ,{
		coefficients[8], coefficients[9], coefficients[10]}};
		double[] B = new double[] {-coefficients[3], -coefficients[7], -1};
		return solveVectorUsingSVD(MatrixUtils.createRealMatrix(A), MatrixUtils.createRealVector(B)).toArray();
	}

	public static RealVector solveVectorUsingSVD(RealMatrix A, RealVector y) {
        // Perform Singular Value Decomposition (SVD)
        SingularValueDecomposition svd = new SingularValueDecomposition(A);
        return svd.getSolver().solve(y);
    }
}
