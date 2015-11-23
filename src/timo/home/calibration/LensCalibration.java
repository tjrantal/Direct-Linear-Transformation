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
	3D DLT (direct linear transformation).
	Depends on Jama library http://math.nist.gov/javanumerics/jama/
	Written by Timo Rantalainen tjrantal@gmail.com
*/

package timo.home.calibration;
import Jama.*;
import java.util.Vector;

public class LensCalibration{
	Matrix calibrationObjectGlobalCoordinates;	//Calibration object global coordinates
	Vector<Matrix> KRt;	//Vector for K, R, and t after calibration
	int u0;	//Optical centre x-coordinate
	int v0;	//Optical centre y-coordinate
	
	/*
		Constructor
		@param calibrationObjectGlobalCoordinates, Nx3 matrix of calibration object global coordinates
		@param u0, optical centre x-coordinate
		@param v0, optical centre y-coordinate		
	*/	
	public LensCalibration(double[][] calibrationObjectGlobalCoordinates, int u0, int v0){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
		this.u0 = u0;
		this.v0 = v0;
		KRt = new resetKRt();
	}

	private Vector<Matrix> resetKRt(){
		temp = new Vector<Matrix>();
		for (int i = 0;i<3;++i){
			temp.add(new Matrix());		
		} 
		return temp;
	}	
	
	/*
		@param digitizedCalibrationObjectCoordinates, Nx2 matrix of digitized calibration object coordinates. insert null for missing coordinate (TO BE IMPLEMENTED)
		@return KRt Vector<Matrix> of camera intrinsic matrix, rotation matrix, and translation vector (in matrix)
	*/
	public Vector<Matrix> calibrate(double[][] digitizedCalibrationObjectCoordinates){
		//IMPLEMENT calibration here
		KRt = new resetKRt(); //Reset KRt
		
		//JATKA TASTA!!
		
		return getCalibration();
	}
	
	public Vector<Matrix> getCalibration(){
		return KRt;
	}

	public void setCalibrationObjectGlobalCoordinates(double[][] calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
	}
	
	public void setCalibrationObjectGlobalCoordinates(Matrix calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = calibrationObjectGlobalCoordinates;
	}
	
		public Vector<Matrix> getDltCoefficients(double[][][] L){
		return getDltCoefficients(this.calibrationObjectGlobalCoordinates,L);
	}
	

	public Vector<Matrix> getDltCoefficients(double[][] calibrationObjectGlobalCoordinates,double[][][] L){
		Vector<Matrix> coefficients = new Vector<Matrix>();
		/*Go through cameras to calibrate*/
		for (int c = 0;c<L.length;++c){
			double[][] B = new double[2*calibrationObjectGlobalCoordinates.length][11];//Matrix for solving DLT-parameters
			double[] C = new double[2*L[c].length]; //Digitized calibrationObject coordinates
			int monta = 0;
			for (int j =0;j<L[c].length;j++){
				for (int i =0;i<2;i++){
					C[monta] = L[c][j][i];
					++monta;
				}
			}
		
			for (int i=0;i<calibrationObjectGlobalCoordinates.length;i++){
				B[2*i][0]	=calibrationObjectGlobalCoordinates[i][0];
				B[2*i][1]	=calibrationObjectGlobalCoordinates[i][1];
				B[2*i][2]	=calibrationObjectGlobalCoordinates[i][2];
				B[2*i][3]	=1;
				B[2*i][4]	=0;
				B[2*i][5]	=0;
				B[2*i][6]	=0;
				B[2*i][7]	=0;
				B[2*i][8]	=-calibrationObjectGlobalCoordinates[i][0]*L[c][i][0];
				B[2*i][9]	=-calibrationObjectGlobalCoordinates[i][1]*L[c][i][0];
				B[2*i][10]	=-calibrationObjectGlobalCoordinates[i][2]*L[c][i][0];
				B[2*i+1][0]	=0;
				B[2*i+1][1]	=0;
				B[2*i+1][2]	=0;
				B[2*i+1][3]	=0;
				B[2*i+1][4]	=calibrationObjectGlobalCoordinates[i][0];
				B[2*i+1][5]	=calibrationObjectGlobalCoordinates[i][1];
				B[2*i+1][6]	=calibrationObjectGlobalCoordinates[i][2];
				B[2*i+1][7]	=1;
				B[2*i+1][8]	=-calibrationObjectGlobalCoordinates[i][0]*L[c][i][1];
				B[2*i+1][9]	=-calibrationObjectGlobalCoordinates[i][1]*L[c][i][1];
				B[2*i+1][10]=-calibrationObjectGlobalCoordinates[i][2]*L[c][i][1];
			}
			//Solve the coefficients
			Matrix A = new Matrix(B);
			Matrix b = new Matrix(C,C.length);
			Matrix coeffs = A.solve(b);
			coefficients.add(coeffs);
		}
		return coefficients;	
	}
	
	
	public Matrix scaleCoordinates(double[][] coordinates){
		return scaleCoordinates(this.dltCoefficients,coordinates);
	}
	
	public Matrix scaleCoordinates(Vector<Matrix> coefficients, double[][] coordinates){
		double[][]	L1 = new double[2*coefficients.size()][3];
		double[]		L2 = new double[2*coordinates.length];
		
		for (int i =0;i<coordinates.length;++i){
			L2[2*i] = coefficients.get(i).get(3,0)- coordinates[i][0];
			L2[2*i+1] = coefficients.get(i).get(7,0)- coordinates[i][1];
		}
		
		for (int i = 0;i<coefficients.size();++i){
			L1[2*i][0]	=coefficients.get(i).get(8,0)*coordinates[i][0]-coefficients.get(i).get(0,0);
			L1[2*i][1]	=coefficients.get(i).get(9,0)*coordinates[i][0]-coefficients.get(i).get(1,0);
			L1[2*i][2]	=coefficients.get(i).get(10,0)*coordinates[i][0]-coefficients.get(i).get(2,0);
			L1[2*i+1][0]	=coefficients.get(i).get(8,0)*coordinates[i][1]-coefficients.get(i).get(4,0);
			L1[2*i+1][1]	=coefficients.get(i).get(9,0)*coordinates[i][1]-coefficients.get(i).get(5,0);
			L1[2*i+1][2]	=coefficients.get(i).get(10,0)*coordinates[i][1]-coefficients.get(i).get(6,0);
		}
		Matrix l1 = new Matrix(L1);
		Matrix l2 = new Matrix(L2,L2.length);
		Matrix result= l1.solve(l2);
		return result;
	}
	
	public Vector<Matrix> getCurrentDltCoefficients(){
		return dltCoefficients;
	}
}
