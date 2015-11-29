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

public class RefineParameters{
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
		double theta = Math.acos((trace(R.getArray())-1d)/2d);
		double[] w;
		if (theta<Math.ulp(0d)){
			w = new double[3];
		}else{
			w = new double[]{
				theta/(2*Math.sin(theta))*(R.get(2,1)-R.get(1,2)),
				theta/(2*Math.sin(theta))*(R.get(0,2)-R.get(2,0)),
				theta/(2*Math.sin(theta))*(R.get(1,0)-R.get(0,1)),
			};
		}
		int noPnts = digit.length;
		int noParam = 11+d.length;
		double[] param = new double[noParam];
		double[][] K_lm = new double[3][3];
		Matrix R_lm;
		double[] w_lm = new double[3];
		double[] t_lm = new double[3];
		double[] d_lm = new double[d.length];
		double[][] J = new double[2*noPnts][noParam];
		double[] dist_lm = new double[2*noPnts];
		double[] dist = new double[2*noPnts];
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
			R_lm = new Matrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
			if (theta>=Math.ulp(0d)){				
				Matrix wh_sk = new Matrix(new double[][]{{0,-w_lm[2],w_lm[1]},{w_lm[2], 0,-w_lm[0]},{-w_lm[1],w_lm[0],0}});
				wh_sk.timesEquals(1/theta);
				//3x3 Rotation matrix
				w = new double[]{
					theta/(2*Math.sin(theta))*(R.get(2,1)-R.get(1,2)),
					theta/(2*Math.sin(theta))*(R.get(0,2)-R.get(2,0)),
					theta/(2*Math.sin(theta))*(R.get(1,0)-R.get(0,1)),
				};
				R_lm.plusEquals(wh_sk.times(Math.sin(theta)));
				R_lm.plusEquals((wh_sk.times(wh_sk)).times(1-Math.cos(theta)));
			}
			//translation vector
			t_lm[0] = t.get(0,0)+delta[8];
			t_lm[1] = t.get(1,0)+delta[9];
			t_lm[2] = t.get(2,0)+delta[10];
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
				XXw[i][calib[i].length] = 1;
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
				xd[i] = xu[i]*(1 + k1*r[i]*r[i] + k2*Math.pow(r[i],4d) + k3*Math.pow(r[i],6d)) + 2*p1*xu[i]*yu[i] + p2*(r[i]*r[i] + 2*xu[i]*xu[i]);
				yd[i] = yu[i]*(1 + k1*r[i]*r[i] + k2*Math.pow(r[i],4d) + k3*Math.pow(r[i],6d)) + p1*(r[i]*r[i] + 2*yu[i]*yu[i]) + 2*p2*xu[i]*yu[i];
				u[i] = K_lm[0][0]*xd[i] + K_lm[0][1]*yd[i] + K_lm[0][2];
			   v[i] = K_lm[1][1]*yd[i] + K_lm[1][2];
			   dist_lm[2*i]	= u[i]-digit[i][0];
			   dist_lm[2*i+1]	= v[i]-digit[i][1];
			}
			//Reprojection error
			rperr_lm = Math.sqrt(dot(dist_lm,dist_lm)/noPnts/2);
			
			if (rperr_lm <= rperr){
				param[0] = K.get(0,0);
				param[1] = K.get(1,1);
				param[2] = K.get(0,2);
				param[3] = K.get(1,2);
				param[4] = K.get(0,1);
				for (int i =0;i<3;++i){
					param[5+i]=w[i];
					param[8+i]=t.get(i,0);
				}
				for (int i =0;i<d.length;++i){
					param[11+i] = d[i];
				}
				K = new Matrix(K_lm);
				R = R_lm.copy();
				t = new Matrix(t_lm,t_lm.length);
				d = Arrays.copyOf(d_lm,d_lm.length);
				rperr = rperr_lm;
				if (n>0 && Math.sqrt(dot(delta,delta)/dot(param,param)) < Math.ulp(0)){
					//Optimization done, exit loop
					break;
				}
				//Update
         	w = Arrays.copyOf(w_lm,w_lm.length);
         	dist  = Arrays.copyOf(dist_lm,dist_lm.length);
         	//Compute the Jacobian
         	for (int i =0;i<noPnts;++i){
         		double[] xxu = new double[]{xu[i],yu[i]};
         		double[] xxd = new double[]{xd[i],yd[i]};
         		//The derivative of a undistorted normalized point
         		Matrix[] dxxu = compute_dxxu(XXw[i], R, t, w);
         		//The derivative of a distorted normalized point
            	Matrix[] dxxd = compute_dxxd(xxu, d, dxxu[0], dxxu[1]);
         		
         		// The derivative of a distotred 2D pixel points
            	Matrix[] dxx = Compute_dxx(xxd, K, d, dxxd[0], dxxd[1], dxxd[2]);
            	//Add the data to the Jacobian
         	}
         	//JATKA TASTA, apufunktiot implementoitu
			}
			
			
		}
			
		//calibration refined
		KdRt = new ArrayList<Matrix>();
		KdRt.add(K);
		KdRt.add(new Matrix(d,d.length));
		KdRt.add(R);
		KdRt.add(t);
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

	/*Return the KdRt ArrayList*/	
	public ArrayList<Matrix> getCalibration(){
		return KdRt;
	}
	
	//Helper functions
	Matrix[] compute_dxxu(double[] XXw,Matrix R,Matrix t,double[] w){
		Matrix XXc = (R.times(new Matrix(XXw,XXw.length).transpose())).plus(t);
		
		Matrix dxxu_dXXc = compute_dxxu_dXXc(XXc);
		Matrix dXXc_dr   = compute_dXXc_dr(XXw);
		Matrix dr_dw     = compute_dr_dw(w);
		Matrix dXXc_dt   = new Matrix(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
		 
		Matrix dxxu_dw = (dxxu_dXXc.times(dXXc_dr)).times(dr_dw);
		Matrix dxxu_dt = dxxu_dXXc.times(dXXc_dt);
		return new Matrix[]{dxxu_dw,dxxu_dt};
	}
	
	Matrix compute_dxxu_dXXc(Matrix XXc){
		double[][] dxxu_dXXc = new double[2][3];
		dxxu_dXXc[0][0] = 1/XXc.get(2,0);  
		dxxu_dXXc[1][0] = 0; 
		dxxu_dXXc[0][1] = 0;
		dxxu_dXXc[1][1] = 1/XXc.get(2,0);
		dxxu_dXXc[0][2] = -XXc.get(0,0)/Math.pow(XXc.get(2,0),2d);
		dxxu_dXXc[1][2] = -XXc.get(1,0)/Math.pow(XXc.get(2,0),2d);;
		return new Matrix(dxxu_dXXc);
	}
	
	Matrix compute_dXXc_dr(double[] XXw){
		double[][] dXXc_dr = new double[3][9];
		for (int i =0;i<XXw.length;++i){
			dXXc_dr[0][i*3] = XXw[i];
			dXXc_dr[1][i*3+1] = XXw[i];
			dXXc_dr[2][i*3+2] = XXw[i];
		}
		return new Matrix(dXXc_dr);
	}
	
	Matrix compute_dr_dw(double[] w){
		double wx = w[0];
		double wy = w[1];
		double wz = w[2];
		double theta = Math.sqrt(dot(w,w));
		double[][] dr_dw = null;
		if (theta < Math.ulp(0)){
			 dr_dw = new double[][]{{ 0 , 0 , 0},
				        {0 , 0 , 1},
				        {0, -1,  0},
				        {0 , 0 ,-1},
				        {0 , 0 , 0},
				        {1 , 0,  0},
				        {0 , 1,  0},
				       {-1 , 0 , 0},
				        {0 , 0 , 0}};
		}else{
			 double wxh = wx/theta;
			 double wyh = wy/theta;
			 double wzh = wz/theta; 

			 double[] dsth_dw = new double[]{wx*Math.cos(theta)/theta, wy*Math.cos(theta)/theta, wz*Math.cos(theta)/theta};
			 double[] domcth_dw = new double[]{wx*Math.sin(theta)/theta, wy*Math.sin(theta)/theta, wz*Math.sin(theta)/theta};

			 double[] dwxh_dw = new double[]{1/theta-Math.pow(wx,2d)/Math.pow(theta,3d), -wx*wy/Math.pow(theta,3d), -wx*wz/Math.pow(theta,3d)};
			 double[] dwyh_dw = new double[]{-wx*wy/Math.pow(theta,3d), 1/theta-Math.pow(wy,2d)/Math.pow(theta,3d), -wy*wz/Math.pow(theta,3d)};
			 double[] dwzh_dw = new double[]{-wx*wz/Math.pow(theta,3d), -wy*wz/Math.pow(theta,3d), 1/theta-Math.pow(wz,2)/Math.pow(theta,3d)};

			 double[] dwxh2pwyh2_dw = new double[]{2*wx/Math.pow(theta,2d)-2*(Math.pow(wx,3d)+Math.pow(wy,2d)*wx)/Math.pow(theta,4d), 2*wy/Math.pow(theta,2d)-2*(Math.pow(wy,3d)+Math.pow(wx,2d)*wy)/Math.pow(theta,4d), -2*(Math.pow(wx,2d)*wz+Math.pow(wy,2d)*wz)/Math.pow(theta,4d)};
			 double[] dwxh2pwzh2_dw = new double[]{2*wx/Math.pow(theta,2d)-2*(Math.pow(wx,3d)+Math.pow(wz,2d)*wx)/Math.pow(theta,4d), -2*(Math.pow(wx,2d)*wy+Math.pow(wz,2d)*wy)/Math.pow(theta,4d), 2*wz/Math.pow(theta,2d)-2*(Math.pow(wz,3d)+Math.pow(wx,2d)*wz)/Math.pow(theta,4d)};
			 double[] dwyh2pwzh2_dw = new double[]{-2*(Math.pow(wy,2d)*wx+Math.pow(wz,2d)*wx)/Math.pow(theta,4d), 2*wy/Math.pow(theta,2d)-2*(Math.pow(wy,3d)+Math.pow(wz,2d)*wy)/Math.pow(theta,4d), 2*wz/Math.pow(theta,2d)-2*(Math.pow(wz,3d)+Math.pow(wy,2d)*wz)/Math.pow(theta,4d)};

			 double wxh2pwyh2 = Math.pow(wxh,2d) + Math.pow(wyh,2d);
			 double wxh2pwzh2 = Math.pow(wxh,2d) + Math.pow(wzh,2d);
			 double wyh2pwzh2 = Math.pow(wyh,2d) + Math.pow(wzh,2d);

			 double omcth = 1-Math.cos(theta);
			 double sth   = Math.sin(theta);

			 dr_dw = new double[][]{{-omcth*dwyh2pwzh2_dw[0] -wyh2pwzh2*domcth_dw[0],-omcth*dwyh2pwzh2_dw[1] -wyh2pwzh2*domcth_dw[1],-omcth*dwyh2pwzh2_dw[2] -wyh2pwzh2*domcth_dw[2]},
				        {sth*dwzh_dw[0] + wzh*dsth_dw[0] + wyh*omcth*dwxh_dw[0] + wxh*omcth*dwyh_dw[0] + wxh*wyh*domcth_dw[0],
				        sth*dwzh_dw[1] + wzh*dsth_dw[1] + wyh*omcth*dwxh_dw[1] + wxh*omcth*dwyh_dw[1] + wxh*wyh*domcth_dw[1],
				        sth*dwzh_dw[2] + wzh*dsth_dw[2] + wyh*omcth*dwxh_dw[2] + wxh*omcth*dwyh_dw[2] + wxh*wyh*domcth_dw[2]},
				       {-sth*dwyh_dw[0] - wyh*dsth_dw[0] + wzh*omcth*dwxh_dw[0] + wxh*omcth*dwzh_dw[0] + wxh*wzh*domcth_dw[0],
				       -sth*dwyh_dw[1] - wyh*dsth_dw[1] + wzh*omcth*dwxh_dw[1] + wxh*omcth*dwzh_dw[1] + wxh*wzh*domcth_dw[1],
				       -sth*dwyh_dw[2] - wyh*dsth_dw[2] + wzh*omcth*dwxh_dw[2] + wxh*omcth*dwzh_dw[2] + wxh*wzh*domcth_dw[2]
				       },
				       {-sth*dwzh_dw[0] - wzh*dsth_dw[0] + wyh*omcth*dwxh_dw[0] + wxh*omcth*dwyh_dw[0] + wxh*wyh*domcth_dw[0],
				       -sth*dwzh_dw[1] - wzh*dsth_dw[1] + wyh*omcth*dwxh_dw[1] + wxh*omcth*dwyh_dw[1] + wxh*wyh*domcth_dw[1], 
				       -sth*dwzh_dw[2] - wzh*dsth_dw[2] + wyh*omcth*dwxh_dw[2] + wxh*omcth*dwyh_dw[2] + wxh*wyh*domcth_dw[2]			       },
				       
				       {-omcth*dwxh2pwzh2_dw[0] - wxh2pwzh2*domcth_dw[0],
				       -omcth*dwxh2pwzh2_dw[1] - wxh2pwzh2*domcth_dw[1],
				       -omcth*dwxh2pwzh2_dw[2] - wxh2pwzh2*domcth_dw[2]},
						{sth*dwxh_dw[0] + wxh*dsth_dw[0] + wzh*omcth*dwyh_dw[0] + wyh*omcth*dwzh_dw[0] + wyh*wzh*domcth_dw[0],
						sth*dwxh_dw[1] + wxh*dsth_dw[1] + wzh*omcth*dwyh_dw[1] + wyh*omcth*dwzh_dw[1] + wyh*wzh*domcth_dw[1],
						sth*dwxh_dw[2] + wxh*dsth_dw[2] + wzh*omcth*dwyh_dw[2] + wyh*omcth*dwzh_dw[2] + wyh*wzh*domcth_dw[2]
						},
				        
				        {sth*dwyh_dw[0] + wyh*dsth_dw[0] + wzh*omcth*dwxh_dw[0] + wxh*omcth*dwzh_dw[0] + wxh*wzh*domcth_dw[0],
				        sth*dwyh_dw[1] + wyh*dsth_dw[1] + wzh*omcth*dwxh_dw[1] + wxh*omcth*dwzh_dw[1] + wxh*wzh*domcth_dw[1],
				        sth*dwyh_dw[2] + wyh*dsth_dw[2] + wzh*omcth*dwxh_dw[2] + wxh*omcth*dwzh_dw[2] + wxh*wzh*domcth_dw[2]},
				        
				       {-sth*dwxh_dw[0] - wxh*dsth_dw[0] + wzh*omcth*dwyh_dw[0] + wyh*omcth*dwzh_dw[0] + wyh*wzh*domcth_dw[0],
				       -sth*dwxh_dw[1] - wxh*dsth_dw[1] + wzh*omcth*dwyh_dw[1] + wyh*omcth*dwzh_dw[1] + wyh*wzh*domcth_dw[1],
				       -sth*dwxh_dw[2] - wxh*dsth_dw[2] + wzh*omcth*dwyh_dw[2] + wyh*omcth*dwzh_dw[2] + wyh*wzh*domcth_dw[2]},
				       {-omcth*dwxh2pwyh2_dw[0] - wxh2pwyh2*domcth_dw[0],
				       -omcth*dwxh2pwyh2_dw[1] - wxh2pwyh2*domcth_dw[1],
				       -omcth*dwxh2pwyh2_dw[2] - wxh2pwyh2*domcth_dw[2]}
				       };
		}
		return new Matrix(dr_dw);
	}
	
	Matrix[] compute_dxxd(double[] xxu,double[] d,Matrix dxxu_dw, Matrix dxxu_dt){ 
		double xu = xxu[0];
		double yu = xxu[1];
		double r  = Math.sqrt(xu*xu + yu*yu);

		double[] dxu_dw = new double[(dxxu_dw.getArray())[0].length];
		double[] dyu_dw = new double[(dxxu_dw.getArray())[1].length];
		for (int i = 0; i<dxu_dw.length;++i){
			dxu_dw[i] = dxxu_dw.get(0,i);
			dyu_dw[i] = dxxu_dw.get(1,i);
		}

		double[] dxu_dt = new double[(dxxu_dt.getArray())[0].length];
		double[] dyu_dt = new double[(dxxu_dt.getArray())[1].length];
		for (int i = 0; i<dxu_dt.length;++i){
			dxu_dt[i] = dxxu_dt.get(0,i);
			dyu_dt[i] = dxxu_dt.get(1,i);
		}

		// The derivative of a radial distortion function
		//   c = 1 + k1*r^2 + k2*r^4 + k3*r^6;
		double c=0;
		double[] dr2_dw = new double[dxu_dw.length];
		double[] dr2_dt = new double[dxu_dt.length];
		double[] dc_dw = new double[dxu_dw.length];
		double[] dc_dt = new double[dxu_dt.length];
		double[] dc_dk = new double[1];
		
		double[] dxr_dw = new double[dxu_dw.length];
		double[] dyr_dw = new double[dxu_dw.length];
		double[] dxr_dt = new double[dxu_dt.length];
		double[] dyr_dt = new double[dxu_dt.length];
		double[] dxr_dk = new double[1];
		double[] dyr_dk = new double[1];

		Matrix dxt_dw = null;
	  	Matrix dyt_dw = null;
		Matrix dxt_dt = null;
		Matrix dyt_dt = null;

		double[] dxt_dp = new double[2];
		double[] dyt_dp = new double[2]; 

		if (d.length >= 1){
			 //First radial distortion coefficients, k1
			 double k1 = d[0];
			 c = 1 + k1*r*r;

			 for (int i = 0;i<dxu_dw.length;++i){
			 	dr2_dw[i] =2*xu*dxu_dw[i] + 2*yu*dyu_dw[i];
			 	dc_dw[i] = k1*dr2_dw[i];
			 }
			 	
			 for (int i = 0;i<dxu_dt.length;++i){
			 	dr2_dt[i] =2*xu*dxu_dt[i] + 2*yu*dyu_dt[i];
			 	dc_dt[i] = k1*dr2_dt[i];
			 }		 
			 dc_dk = new double[]{r*r};
		}

		 if (d.length >= 2){
			  //Second radial distortion coefficients, k2
			  double k2 = d[1];
			  c = c + k2*Math.pow(r,4d);
			  for (int i = 0; i<dc_dw.length;++i){
			  		dc_dw[i] = dc_dw[i]+k2*2*r*r*dr2_dw[i];
			  }
			  for (int i = 0; i<dc_dt.length;++i){
			  		dc_dt[i] = dc_dt[i] + k2*2*r*r*dr2_dt[i];
			  }
			  dc_dk = new double[]{dc_dk[0], Math.pow(r,4d)};
		 }

		 if (d.length == 5){
			  //Third radial distortion coefficients, k3
			  double k3 = d[4];
			  c = c + k3*Math.pow(r,6d);
			  for (int i = 0; i<dc_dw.length;++i){
			  		dc_dw[i] = dc_dw[i]+k3*3*Math.pow(r,4d)*dr2_dw[i];
			  }
			  for (int i = 0; i<dc_dt.length;++i){
			  		dc_dt[i] = dc_dt[i] + k3*3*Math.pow(r,4d)*dr2_dt[i];
			  }
			  dc_dk = new double[]{dc_dk[0],dc_dk[1], Math.pow(r,6d)};
		 }

		 if (d.length > 0){
			  //The derivative of a radially distorted normalized point
			  for (int i =0;i<dc_dw.length;++i){
				  dxr_dw[i] = c*dxu_dw[i] + xu*dc_dw[i];
				  dyr_dw[i] = c*dyu_dw[i] + yu*dc_dw[i];
				}
				for (int i =0;i<dc_dt.length;++i){
				  	dxr_dt[i] = c*dxu_dt[i] + xu*dc_dt[i];
				  	dyr_dt[i] = c*dyu_dt[i] + yu*dc_dt[i];
				}
				dxr_dk = new double[dc_dk.length];
				dyr_dk = new double[dc_dk.length];
				for (int i =0;i<dc_dk.length;++i){
				  	dxr_dk[i] = xu*dc_dk[i];
				  	dyr_dk[i] = yu*dc_dk[i];
			  	}
		 }


		 //The derivative of a tangentially distorted normalized point
		 if (d.length >= 4){
			  double p1 = d[2];
			  double p2 = d[3];
			  
			  Matrix dxt_dxxu = new Matrix(new double[]{2*p1*yu+6*p2*xu, 2*p1*xu+2*p2*yu},2).transpose();
			  Matrix dyt_dxxu = new Matrix(new double[]{2*p1*xu+2*p2*yu, 6*p1*yu+2*p2*xu},2).transpose();
				
				double[][] tw = new double[dxu_dw.length][2];
				double[][] tt = new double[dxu_dt.length][2];	
				for (int i =0;i<tw.length;++i){
					tw[i][0] = dxu_dw[i];
					tw[i][1] = dyu_dw[i];
				}
				for (int i =0;i<tt.length;++i){
					tt[i][0] = dxu_dt[i];
					tt[i][1] = dyu_dt[i];
				}
			  dxxu_dw = new Matrix(tw);
			  dxxu_dt = new Matrix(tt);

			  dxt_dw = dxt_dxxu.times(dxxu_dw);
			  dyt_dw = dyt_dxxu.times(dxxu_dw);

			  dxt_dt = dxt_dxxu.times(dxxu_dt);
			  dyt_dt = dyt_dxxu.times(dxxu_dt);

			  dxt_dp = new double[]{2*xu*yu, r*r+2*xu*xu};
			  dyt_dp = new double[]{r*r+2*yu*yu, 2*xu*yu};
		 }


		//The derivative of a distorted normalized point
		double[]	dxd_dw = dxr_dw;
		 double[] dyd_dw = dyr_dw;

		 double[] dxd_dt = dxr_dt;
		 double[] dyd_dt = dyr_dt;

		 double[] dxd_dd = dxr_dk;
		 double[] dyd_dd = dyr_dk;

		if (d.length > 2){
			for (int i = 0;i<dxd_dw.length;++i){
			 dxd_dw[i] += dxt_dw.get(i,0);
			 dyd_dw[i] += dyt_dw.get(i,0);
			}
			for (int i = 0;i<dxd_dt.length;++i){
			 dxd_dt[i] += dxt_dt.get(i,0);
			 dyd_dt[i] += dyt_dt.get(i,0);
			}

			dxd_dd = new double[]{dxr_dk[0], dxr_dk[1], dxt_dp[0], dxt_dp[1], dxr_dk[2]};
			dyd_dd = new double[]{dyr_dk[0], dyr_dk[1], dyt_dp[0], dyt_dp[1], dyr_dk[2]};
		}

		double[][] dxxd_dw = new double[2][dxd_dw.length];
		for (int i = 0;i<dxd_dw.length;++i){
			dxxd_dw[0][i] = dxd_dw[i];
			dxxd_dw[1][i] = dyd_dw[i];
		}
		
		double[][] dxxd_dt = new double[2][dxd_dt.length];
		for (int i = 0;i<dxd_dw.length;++i){
			dxxd_dt[0][i] = dxd_dt[i];
			dxxd_dt[1][i] = dyd_dt[i];
		}
		
		double[][] dxxd_dd = new double[2][dxd_dd.length];
		for (int i = 0;i<dxd_dw.length;++i){
			dxxd_dd[0][i] = dxd_dd[i];
			dxxd_dd[1][i] = dyd_dd[i];
		}
		return new Matrix[]{new Matrix(dxxd_dw),new Matrix(dxxd_dt),new Matrix(dxxd_dd)};
	}
	
	Matrix[] Compute_dxx(double[] xxd, Matrix K, double[] d, Matrix dxxd_dw,Matrix dxxd_dt,Matrix dxxd_dd){
		double xd = xxd[0];
		double yd = xxd[1];

		double fx = K.get(0,0);
		double fy = K.get(1,1);
		double u0 = K.get(0,2);
		double v0 = K.get(1,2);
		double s  = K.get(0,1);

		double[] dxd_dw = (dxxd_dw.getArray())[0];
		double[] dyd_dw = (dxxd_dw.getArray())[1];

		double[] dxd_dt = (dxxd_dt.getArray())[0];
		double[] dyd_dt = (dxxd_dt.getArray())[1];

		double du_dfx = xd;
		double du_dfy =  0;
		double du_du0 =  1;
		double du_dv0 =  0;
		double dv_dfx =  0;
		double dv_dfy = yd;
		double dv_du0 =  0;
		double dv_dv0 =  1; 
		double du_ds,dv_ds;
		if (s == 0){
			 du_ds = 0;
			 dv_ds = 0;
		}else{
			 du_ds = yd;
			 dv_ds =  0;
		}
		
		double[] dxd_dd = (dxxd_dd.getArray())[0];
		double[] dyd_dd = (dxxd_dd.getArray())[1];
		double[] du_dd = new double[dxd_dd.length];
		double[] dv_dd = new double[dxd_dd.length];
		for (int i = 0;i<dxd_dd.length;++i){
			 du_dd[i] = fx*dxd_dd[i] + s*dyd_dd[i];
			 dv_dd[i] = fy*dyd_dd[i];
		 }


		double[] du_dk = new double[]{du_dfx, du_dfy, du_du0, du_dv0, du_ds};
		double[] dv_dk = new double[]{dv_dfx, dv_dfy, dv_du0, dv_dv0, dv_ds};

		double[] du_dw = new double[dxd_dw.length];
		double[] dv_dw = new double[dxd_dw.length];
		for (int i = 0;i<dxd_dd.length;++i){
			 du_dw[i] = fx*dxd_dw[i] + s*dyd_dw[i];
			 dv_dw[i] = fy*dyd_dw[i];
		 }

		double[] du_dt = new double[dxd_dt.length];
		double[] dv_dt = new double[dxd_dt.length];
		for (int i = 0;i<dxd_dd.length;++i){
			 du_dt[i] = fx*dxd_dt[i] + s*dyd_dt[i];
			 dv_dt[i] = fy*dyd_dt[i];
		 }

		double[][] dxx_dk = new double[][]{du_dk, dv_dk};
		double[][] dxx_dw = new double[][]{du_dw, dv_dw};
		double[][] dxx_dt = new double[][]{du_dt, dv_dt};
		double[][] dxx_dd = new double[][]{du_dd, dv_dd};
		return new Matrix[]{new Matrix(dxx_dk),new Matrix(dxx_dw),new Matrix(dxx_dt),new Matrix(dxx_dd)};
	}
}
