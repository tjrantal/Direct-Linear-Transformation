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

package timo.home.refine;
import Jama.*;
import java.util.*;

public class LensCalibration{
	Matrix calibrationObjectGlobalCoordinates;	//Calibration object global coordinates
	ArrayList<Matrix> KdRt;	//ArrayList for K, R, and t after calibration
	double[][] calib;	//Calibration object global coordinates
	double[][] digit;	//Digitized object coordinates
	Matrix K;			//Camera matrix
	double[] d;			//distortion coefficients
	Matrix R;			//Rotation matrix
	Matrix t;			//Translation vector
	
	/*
		Constructor
		@param calib, Nx3 matrix of calibration object global coordinates
		@param digit, Nx2 matrix of calibration object digitized coordinates
		@param K, Camera matrix from LensCalibration
		@param d, 2x1, or 5x1 array of distortion coefficients (can be array of zeroes)
		@param R, Rotation matrix from LensCalibration
		@param t, Translation matrix from LensCalibration		
	*/	
	public RefineParameters(double[][] calib,double[][] digit, Matrix K, double[] d, Matrix R, Matrix t){
		this.calib = calib;
		this.digit = digit;
		this.K = K;
		this.d = d;
		this.R = R;
		this.t = t;
		KdRt = null;
		refine();
	}


	/*
		Refines calibration with LM(Levenberg-Marquardt) nonlinear least squares algorithm. Ported from http://www.daesik80.com/matlabfns/matlabfns.htm
		@return KdRt ArrayList<Matrix> of camera intrinsic matrix, lens distortion coefficients, rotation matrix, and translation vector (in matrix)
	*/
	public ArrayList<Matrix> refine(){
		//IMPLEMENT calibration here
		double centred[][] = new double[digitizedCalibrationObjectCoordinates.length][2];
		
		for (int i = 0; i<centred.length;++i){
			centred[i][0] = digitizedCalibrationObjectCoordinates[i][0]-u0;
			centred[i][1] = digitizedCalibrationObjectCoordinates[i][1]-v0;
		}
		/*Prep the matrix for SVD*/
		double[][] M = new double[centred.length][8];
		double[][] calib = calibrationObjectGlobalCoordinates.getArray();
		for (int i = 0; i<centred.length;++i){
			M[i][0] = centred[i][0]*calib[i][0];
			M[i][1] = centred[i][0]*calib[i][1];
			M[i][2] = centred[i][0]*calib[i][2];
			M[i][3] = centred[i][0];
			M[i][4] = -centred[i][1]*calib[i][0];
			M[i][5] = -centred[i][1]*calib[i][1];
			M[i][6] = -centred[i][1]*calib[i][2];
			M[i][7] = -centred[i][1];
		}
		SingularValueDecomposition svd = new SingularValueDecomposition(new Matrix(M));
		Matrix V = svd.getV();
		double scale = 1/Math.sqrt(Math.pow(V.get(0,7),2d)+Math.pow(V.get(1,7),2d)+Math.pow(V.get(2,7),2d));
		double alpha = scale*Math.sqrt(Math.pow(V.get(4,7),2d)+Math.pow(V.get(5,7),2d)+Math.pow(V.get(6,7),2d));
		
		if (scale*centred[0][1]*(V.get(0,7)*calib[0][0] + V.get(1,7)*calib[0][1] + V.get(2,7)*calib[0][2] + V.get(3,7)) < 0){
    		scale = -scale;
    	}
		
		double[][] r = new double[3][3];
		r[0][0] = scale*V.get(4,7)/alpha;
		r[0][1] = scale*V.get(5,7)/alpha;
		r[0][2] = scale*V.get(6,7)/alpha;
		r[1][0] = scale*V.get(0,7);
		r[1][1] = scale*V.get(1,7);
		r[1][2] = scale*V.get(2,7);
		double[] t = new double[3];
		t[0] = scale*V.get(7,7)/alpha;
		t[1] = scale*V.get(3,7);
		r = getCross(r);
		SingularValueDecomposition svd2 = new SingularValueDecomposition(new Matrix(r));
		Matrix R = svd2.getU().times(svd2.getV().transpose());
		
		//Compute fx, fy, and tz
		double[] Xc = new double[calib.length];
		double[] Yc = new double[calib.length];
		double[] Zc_tz = new double[calib.length];
		for (int i =0;i<calib.length;++i){
			Xc[i]		= R.get(0,0)*calib[i][0] + R.get(0,1)*calib[i][1] + R.get(0,2)*calib[i][2] + t[0];
			Yc[i] 	=  R.get(1,0)*calib[i][0] + R.get(1,1)*calib[i][1] + R.get(1,2)*calib[i][2] + t[1];
			Zc_tz[i]	= R.get(2,0)*calib[i][0] + R.get(2,1)*calib[i][1] + R.get(2,2)*calib[i][2];
		}
		
		double[][] A = new double[centred.length*2][2];
		double[] B = new double[centred.length*2];
		
		for (int i = 0; i<centred.length;++i){
			A[2*i][0]	= Xc[i];
			A[2*i][1]	= -centred[i][0];
			A[2*i+1][0] = Yc[i]/alpha;
			A[2*i+1][1] = -centred[i][1];
			
			B[2*i]		= centred[i][0]*Zc_tz[i];
			B[2*i+1]		= centred[i][1]*Zc_tz[i];
		}
		Matrix Y = new Matrix(A).inverse().times(new Matrix(B,B.length));

		double[][] K = new double[3][3];
		K[0][0] = Y.get(0,0);
		K[0][2] = u0;
		K[1][1] = Y.get(0,0)/alpha;
		K[1][2] = v0;
		K[2][2] = 1;
		t[2] = Y.get(1,0);
		KRt = new ArrayList<Matrix>();
		KRt.add(new Matrix(K));
		KRt.add(R);
		KRt.add(new Matrix(t,t.length));
		return getCalibration();
	}
	
	/*calculate cross-product from the first two rows onto the third*/
	double[][] getCross(double[][] temp){
		temp[2][0] = temp[0][1]*temp[1][2]-temp[1][1]*temp[0][2];
		temp[2][1] = temp[0][2]*temp[1][0]-temp[1][2]*temp[0][0];
		temp[2][2] = temp[0][0]*temp[1][1]-temp[1][0]*temp[0][1];
		return temp;
	}
	
	public ArrayList<Matrix> getCalibration(){
		return KRt;
	}

	public void setCalibrationObjectGlobalCoordinates(double[][] calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = new Matrix(calibrationObjectGlobalCoordinates);
	}
	
	public void setCalibrationObjectGlobalCoordinates(Matrix calibrationObjectGlobalCoordinates){
		this.calibrationObjectGlobalCoordinates = calibrationObjectGlobalCoordinates;
	}
}
