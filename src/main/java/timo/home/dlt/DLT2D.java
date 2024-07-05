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

    Copyright (C) 2012 Timo Rantalainen
*/
/*
	2D DLT (direct linear transformation).
	Depends on Jama library http://math.nist.gov/javanumerics/jama/
	Written by Timo Rantalainen tjrantal@gmail.com
*/

package timo.home.dlt;
import Jama.*;

public class DLT2D{
	Matrix calibrationObjectGlobalCoordinates;
	Matrix dltCoefficients;
	public DLT2D(){}
	public DLT2D(double[][] calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
	}
	
	public DLT2D(double[][] calibrationObjectGlobalCoordinates, double[][] digitizedCalibrationObjectCoordinates){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
		this.dltCoefficients = getDltCoefficients(calibrationObjectGlobalCoordinates,digitizedCalibrationObjectCoordinates);		
	}
	
	public void setCalibrationObjectGlobalCoordinates(double[][] calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
	}
	
	public void setCalibrationObjectGlobalCoordinates(Matrix calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = calibrationObjectGlobalCoordinates;
	}
	
		public Matrix getDltCoefficients(double[][] L){
		return getDltCoefficients(this.calibrationObjectGlobalCoordinates,L);
	}
	
	public Matrix getDltCoefficients(Matrix calibrationObjectGlobalCoordinates,double[][] L){
		return getDltCoefficients(calibrationObjectGlobalCoordinates.getArray(),L);
	}
	
	public Matrix getDltCoefficients(double[][] calibrationObjectGlobalCoordinates,double[][] L){
		double[][] B = new double[2*calibrationObjectGlobalCoordinates.length][8];//Matrix for solving DLT-parameters
		double[] C = new double[2*L.length]; //Digitized calibrationObject coordinates
		int monta = 0;
		for (int j =0;j<L.length;j++){
			for (int i =0;i<2;i++){
				C[monta] = L[j][i];
				++monta;
			}
		}

		for (int i= 0;i<calibrationObjectGlobalCoordinates.length;++i){
		  B[2*i][0]			= calibrationObjectGlobalCoordinates[i][0]; 
		  B[2*i][1]			= calibrationObjectGlobalCoordinates[i][1]; 
		  B[2*i][2]  		= 1;
		  B[2*i][6]		 	=-calibrationObjectGlobalCoordinates[i][0]*L[i][0];
		  B[2*i][7]			=-calibrationObjectGlobalCoordinates[i][1]*L[i][0];
		  B[2*i+1][3]		= calibrationObjectGlobalCoordinates[i][0];
		  B[2*i+1][4]		= calibrationObjectGlobalCoordinates[i][1];
		  B[2*i+1][5]		= 1;
		  B[2*i+1][6]		=-calibrationObjectGlobalCoordinates[i][0]*L[i][1];
		  B[2*i+1][7]		=-calibrationObjectGlobalCoordinates[i][1]*L[i][1];
		}

		//Solve the coefficients
		Matrix A = new Matrix(B);
		Matrix b = new Matrix(C,C.length);
		Matrix coefficients = A.solve(b);
		return coefficients;	
	}

	public Matrix scaleCoordinates(double[] coordinates){
		return scaleCoordinates(this.dltCoefficients,coordinates);
	}
	
	public Matrix scaleCoordinates(Matrix coefficients, double[] coordinates){
		double[][]	L1 = new double[2][2];
		double[]		L2 = new double[2];
		L1[0][0]	= coefficients.get(0,0)-coordinates[0]*coefficients.get(6,0);
		L1[0][1]	= coefficients.get(1,0)-coordinates[0]*coefficients.get(7,0);
		L1[1][0]	= coefficients.get(3,0)-coordinates[1]*coefficients.get(6,0);
		L1[1][1]	= coefficients.get(4,0)-coordinates[1]*coefficients.get(7,0);
		L2[0]		= coordinates[0]-coefficients.get(2,0);
		L2[1]		= coordinates[1]-coefficients.get(5,0);

		Matrix l1 = new Matrix(L1);
		Matrix l2 = new Matrix(L2,L2.length);
		Matrix result= l1.solve(l2);
		return result;
	}
	
	public Matrix getCurrentDltCoefficients(){
		return dltCoefficients;
	}
}
