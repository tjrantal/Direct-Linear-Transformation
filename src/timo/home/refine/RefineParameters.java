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
		//Start refining calibration
		double[][] r = R.getArray();
		double theta = Math.acos((trace(r)-1d)/2d);
		double[] w;
		if (theta<Math.ulp(0d)){
			w = new double[3];
		}else{
			w = {
				theta/(2*Math.sin(theta))*(R.get(2,1)-R.get(1,2)),
				theta/(2*Math.sin(theta))*(R.get(0,2)-R.get(2,0)),
				theta/(2*Math.sin(theta))*(R.get(1,0)-R.get(0,1)),
			};
		}
		int noPnts = digit.length
		int noParam = 11+d.length;
		double[] param = new double(noParam);
		double[][] K_lm = new double[3][3];
		Matrix R_lm;
		double[] w_lm = new double[3];
		double[] t_lm = new double[3];
		double[] d_lm = new double[d.length];
		double[][] J = new double[2*noPnts][noParam];
		double[] dist_lm = new double[2*noPnts];
		double[] delta = new double[2*noParam];
		double k1,k2,k3,p1,p2,rperr_lm;
		double rperr = Double.POSITIVE_INFINITY; 
		
		for (int n = 0;n<30;++n){
			//Camera intrinsic parameters
			K_lm[0][0] = K.get(0,0)+delta[0];	//fx
			K_lm[1][1] = K.get(1,1)+delta[1];	//fy
			K_lm[0][2] = K.get(0,2)+delta[2];	//u0
			K_lm[1][2] = K.get(1,2)+delta[3];	//v0
			K_lm[0][1] = K.get(0,1)+delta[4];	//s
			K_lm[2][2] = 1;
			
			//3x1 vector of Rodigrues representation into 3x3 rotation matrix
			w_lm[0] = w[0]+delta[5];
			w_lm[1] = w[1]+delta[6];
			w_lm[2] = w[2]+delta[7];
			theta = Math.sqrt(w_lm[0]*w_lm[0]+w_lm[1]*w_lm[1]+w_lm[2]*w_lm[2]);
			R_lm = new Matrix({{1,0,0},{0,1,0},{0,0,1}});
			if (theta>=Math.ulp(0d)){				
				Matrix wh_sk = new Matrix({{0,-w_lm[2],w_lm[1]},{w_lm[2], 0,-w_lm[0]},{-w_lm[1],w_lm[0],0}});
				wh_sk.timesEquals(1/theta);
				//3x3 Rotation matrix
				w = new double{
					theta/(2*Math.sin(theta))*(R.get(2,1)-R.get(1,2)),
					theta/(2*Math.sin(theta))*(R.get(0,2)-R.get(2,0)),
					theta/(2*Math.sin(theta))*(R.get(1,0)-R.get(0,1)),
				};
				R_lm.plusEquals(wh_sk.times(Math.sin(theta)));
				R_lm.plusEquals((wh_sk.times(wh_sk)).times(1-Math.cos(theta)));
			}
			//translation vector
			t_lm[0] = t[0]+delta[8];
			t_lm[1] = t[1]+delta[9];
			t_lm[2] = t[2]+delta[10];
			//3D points presented with respect to the camera coordinates
			double[][] rotTrans = new double[3][4];
			for (int i = 0; i<3;++i){
				for (int j = 0; j<3;++j){
					rotTrans[i][j] = R_lm.get(i,j);
				}
			}
			for (int i = 0; i<3;++i){
				rotTrans[i][3] = t_lm[i];
			}
			double[][] XXw = new double[calib.length][calib[0].length+1];
			for (int i = 0;i<calib.length;++i){
				for (int j = 0; j<calib[i].length;++j){
					XXw[i][j] = calib[i][j];
				}
				XXw[i][calib.length] = 1;
			}
			Matrix XXc = new Matrix(rotTrans).times(new Matrix(XXw));
			//Various d lengths should be implemented here
			 k1 = d[0] + delta[11];
			 k2 = d[1] + delta[12]; 
          p1 = d[2] + delta[13];
          p2 = d[3] + delta[14]; 
          k3 = d[4] + delta[15];
          d_lm[0] = k1;
          d_lm[1] = k2;
          d_lm[2] = p1;
          d_lm[3] = p2;
          d_lm[4] = k3;	
			double[] xu = new double[XXw.length];
			double[] yu = new double[XXw.length];
			double[] xd = new double[XXw.length];
			double[] yd = new double[XXw.length];
			double[] u = new double[XXw.length];
			double[] v = new double[XXw.length];
			double[] r = new double[xu.length];
			for (int i = 0; i<xu.length;++i){
				//undistorted coordinates
				xu[i] = XXc.get(0,i)/XXc.get(2,i);
				yu[i] = XXc.get(1,i)/XXc.get(2,i);
				r[i] = Math.sqrt(xu[i]*xu[i]+yu[i]*yu[i]);
				//distorted coordinates
				xd[i] = xu[i]*(1 + k1*r[i]*r[i] + k2*Math.pow(r[i],4d) + k3*Math.pow(r[i],6d)) + 2*p1*xu[i]*yu[i] + p2*(r[i]*r[i] + 2*xu[i]*xy[i]);
				yd[i] = yu[i]*(1 + k1*r[i]*r[i] + k2*Math.pow(r[i],4d) + k3*Math.pow(r[i],6d)) + p1*(r[i]*r[i] + 2*yu[i]*yu[i]) + 2*p2*xu[i]*yu[i];
				u[i] = K_lm[0][0]*xd[i] + K_lm[0][1]*yd[i] + K_lm[0][2];
			   v[i] = K_lm[1][1]*yd[i] + K_lm[1][2];
			   dist_lm[2*i]	= u[i]-x[i][0];
			   dist_lm[2*i+1]	= v[i]-x[i][1];
			}
			//Reprojection error
			rperr_lm = Math.sqrt(dot(dist_lm,dist_lm)/noPnts/2);
			
			//JATKA TASTA
			
		}
			
		//calibration refined
		KdRt = new ArrayList<Matrix>();
		KdRt.add(new Matrix(K));
		KdRt.add(new Matrix(d,d.length));
		KdRt.add(R);
		KdRt.add(new Matrix(t,t.length));
		return getCalibration();
	}
	
	/*Get trace of a N x N matrix*/
	double trace(double[][] matrix){
		double temp = 0;
		for (int i = 0;i<matrix.length;++i){
			temp+=matrix[i][i];
		}
		return temp;
	}
	
	/*calculate cross-product from the first two rows onto the third*/
	double[][] getCross(double[][] temp){
		temp[2][0] = temp[0][1]*temp[1][2]-temp[1][1]*temp[0][2];
		temp[2][1] = temp[0][2]*temp[1][0]-temp[1][2]*temp[0][0];
		temp[2][2] = temp[0][0]*temp[1][1]-temp[1][0]*temp[0][1];
		return temp;
	}
	
	/*Calculate the dot product between a and b*/
	double dot(double[] a, double[] b){
		double temp = 0;
		for (int i = 0;i<a.length;++i){
			temp+=a[i]*b[i];
		}
		return temp;
	}
	
	public ArrayList<Matrix> getCalibration(){
		return KdRt;
	}

}
