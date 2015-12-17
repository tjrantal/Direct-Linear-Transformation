package timo.home.undristort;

import Jama.*;
import java.awt.image.*;
import timo.home.imagePanel.ImagePanel;

public class Undistort{
	public BufferedImage ubi;	//Undistorted buffered image
	Matrix K;
	Matrix d;
	Matrix invK;
	/*Constructor*/
	public Undistort(ImagePanel ip, Matrix K, Matrix d){
		this.K = K;
		this.d = d;
		invK = K.inverse();
		
		BufferedImage bi = ip.getOrig();
		int width =bi.getWidth();
		int height = bi.getHeight();
		byte[][][] pixArr = getPixelArray(bi);
		byte[][][] undistorted = new byte[height][width][3];
		double x_u,y_u,radius,radial;
		double[] xx_d = new double[3];
		double k1 = d.get(0,0);
		double k2 = d.get(1,0);
		double k3 = d.get(4,0);
		double p1 = d.get(2,0);
		double p2 = d.get(3,0);
		double[][] xxTemp = new double[3][width*height];
		for (int r = 0;r<height;++r){
			for (int c = 0; c<width;++c){
				xxTemp[0][c+r*width] = c;
				xxTemp[1][c+r*width] = r;
				xxTemp[2][c+r*width] = 1d;
			}
		}
		Matrix xx_u = invK.times(new Matrix(xxTemp));
		double[][] xxdTemp = new double[3][width*height];
		double[][] xx_up= xx_u.getArray();
		for (int i = 0;i<width*height;++i){
			x_u = xx_up[0][i];
			y_u = xx_up[1][i];
			radius = Math.sqrt(x_u*x_u + y_u*y_u);
			radial = (1 + k1*radius*radius + k2*radius*radius*radius*radius  + k3*radius*radius*radius*radius*radius*radius);
			xxdTemp[0][i] = radial*x_u + 2*p1*x_u*y_u + p2*(radius*radius + 2*x_u*x_u);
			xxdTemp[1][i] = radial*y_u + p1*(radius*radius + 2*y_u*y_u) + 2*p2*x_u*y_u;
			xxdTemp[2][i] = 1d;
		}
		Matrix xx_dim = K.times(new Matrix(xxdTemp));
		double[][] xx_dimp = xx_dim.getArray();
		for (int r = 0;r<height;++r){
			for (int c = 0; c<width;++c){
				undistorted[r][c] = getBilinearInterpolatedPixel(xx_dimp[0][c+r*width],xx_dimp[1][c+r*width],pixArr);
				//undistorted[r][c] = pixArr[r][c];								
			}
		}
		ubi = ip.createImageFromBytes(undistorted);	
	}
	
	public double[][] projectKnownPoints(double[][] coordinates, Matrix R, Matrix t){
			//Transpose the coordinates from Nx2 to 2xN (and add a row of 1s) 
			double[][] XXw = new double[coordinates[0].length+1][coordinates.length];
			//System.out.println("XXw " +XXw.length+" "+XXw[0].length);
			for (int i = 0;i<coordinates.length;++i){
				for (int j = 0; j<coordinates[i].length;++j){
					XXw[j][i] = coordinates[i][j];
				}
				XXw[coordinates[i].length][i] = 1;
			}
			//3D points presented with respect to the camera coordinates
			double[][] rotTrans = new double[3][4];
			for (int i = 0; i<3;++i){
				for (int j = 0; j<3;++j){
					rotTrans[i][j] = R.get(i,j);
				}
			}
			for (int i = 0; i<3;++i){
				rotTrans[i][3] = t.get(i,0);
			}

			//System.out.println("RotTrans "+(n+1));
			//new Matrix(rotTrans).print(4,4);
			Matrix XXc = new Matrix(rotTrans).times(new Matrix(XXw));
			double[][] cameraCoords = new double[3][coordinates.length];
			for (int i = 0; i<coordinates.length;++i){
				cameraCoords[0][i]= XXc.get(0,i)/XXc.get(2,i);
				cameraCoords[1][i] = XXc.get(1,i)/XXc.get(2,i);
				cameraCoords[2][i]= 1d;
			}
			Matrix temp = K.times(new Matrix(cameraCoords));
			double[][] tempP = temp.getArray();
			double[][] returnCoords = new double[coordinates.length][2];
			for (int i = 0;i<coordinates.length;++i){
					returnCoords[i][0] = tempP[0][i];
					returnCoords[i][1] = tempP[1][i];
			}
			
			return returnCoords;	
	}
	
	/*
		Has to be solved iteratively, used Gauss-Newton method (or could be just Newton's)
		Used Sage to get the partial derivatives
		
		k1,k2,k3,p1,p2,u,v,x,y = var('k1, k2, k3, p1, p2, u, v, x, y')
		0==u*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+2*p1*u*v+p2*((u^2+v^2)+2*u^2)-x
		0==v*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+p1*((u^2+v^2)+2*v^2)+2*p2*u*v-y

		diff(0==u*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+2*p1*u*v+p2*((u^2+v^2)+2*u^2)-x,u)
		diff(0==u*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+2*p1*u*v+p2*((u^2+v^2)+2*u^2)-x,v)
		diff(0==v*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+p1*((u^2+v^2)+2*v^2)+2*p2*u*v-y,u)
		diff(0==v*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+p1*((u^2+v^2)+2*v^2)+2*p2*u*v-y,v)


		
	*/
	public double[][] undistortCoordinates(double[][] coordinatesIn){
		double[][] undistorted = new double[coordinatesIn.length][2];
		/*Exclude missing points*/
		//Check for missing coordinates
		int validCoords = 0;
		for (int i = 0; i<coordinatesIn.length;++i){
			if (coordinatesIn[i] != null){
				++validCoords;
			}else{
				undistorted[i] = null;	//Set invalid output to null
			}
		}
		//Indices for valid coordinates
		int[] validIndices = new int[validCoords];
		validCoords = 0;
		for (int i = 0; i<coordinatesIn.length;++i){
			if (coordinatesIn[i] != null){
				validIndices[validCoords++] = i;
			}	
		}
		
		double[][] coordinates = new double[validIndices.length][2];
		
		for (int i = 0; i<validIndices.length;++i){
			coordinates[i][0] = coordinatesIn[validIndices[i]][0];
			coordinates[i][1] = coordinatesIn[validIndices[i]][1];
		}
		
		
		double[] xx_d = new double[3];
		double k1 = d.get(0,0);
		double k2 = d.get(1,0);
		double k3 = d.get(4,0);
		double p1 = d.get(2,0);
		double p2 = d.get(3,0);
		double[][] xxTemp = new double[3][coordinates.length];
		for (int i = 0;i<coordinates.length;++i){
				xxTemp[0][i] = coordinates[i][0];
				xxTemp[1][i] = coordinates[i][1];
				xxTemp[2][i] = 1d;
		}
		Matrix xx_u = invK.times(new Matrix(xxTemp));
		double[][] xxdTemp = new double[3][coordinates.length];
		double[][] xx_up= xx_u.getArray();
		double x_u,y_u,x_d,y_d,diffXY,x_t,y_t;
		int iterations;
		Matrix xy1;
		for (int i = 0;i<coordinates.length;++i){
			x_u = xx_up[0][i];
			y_u = xx_up[1][i];
			x_d = x_u;
			y_d = y_u;
			diffXY = 100;
			iterations = 0;
			while (diffXY > Math.ulp(1d) && iterations < 100){
				xy1 = new Matrix(new double[]{x_u,y_u},2).minus(new Matrix(jac(x_u,y_u,k1,k2,p1,p2,k3)).inverse().times(new Matrix(funct(x_u,y_u,x_d,y_d,k1,k2,p1,p2,k3),2)));
				diffXY = Math.sqrt(Math.pow(xy1.get(0,0)-x_u,2d)+Math.pow(xy1.get(1,0)-x_u,2d));
				x_u = xy1.get(0,0);
				y_u = xy1.get(1,0);
				++iterations;
			}
			
			xxdTemp[0][i] = x_u;
			xxdTemp[1][i] = y_u;
			xxdTemp[2][i] = 1d;
		}
		Matrix xx_dim = K.times(new Matrix(xxdTemp));
		double[][] xx_dimp = xx_dim.getArray();
		for (int i = 0;i<coordinates.length;++i){
				undistorted[validIndices[i]][0] = xx_dimp[0][i];
				undistorted[validIndices[i]][1] = xx_dimp[1][i];
		}
		return undistorted;
	}
	
	/**func*/
	double[] funct(double u, double v, double x, double y, double k1, double k2, double p1, double p2, double k3){
		double rSq = Math.pow(u,2d)+Math.pow(v,2d);
		double radial = (1+k1*(rSq)+k2*Math.pow((rSq),2d)+k3*Math.pow((rSq),3d));
		return new double[]{
			u*radial+2*p1*u*v+p2*((rSq)+2*Math.pow(u,2d))-x,
		 	v*radial+p1*((rSq)+2*Math.pow(v,2d))+2*p2*u*v-y
		};
	}
	
	/**Jacobian of func*/
	double[][] jac(double u, double v, double k1, double k2, double p1, double p2, double k3){
		double rSq = Math.pow(u,2d)+Math.pow(v,2d);
		return new double[][] { 
			{
				Math.pow(rSq,3d)*k3 + Math.pow(rSq,2d)*k2 + (rSq)*k1 + 2*(3*Math.pow(rSq,2d)*k3*u + 2*(rSq)*k2*u + k1*u)*u + 6*p2*u + 2*p1*v + 1, 
				2*(3*Math.pow(rSq,2d)*k3*v + 2*(rSq)*k2*v + k1*v)*u + 2*p1*u + 2*p2*v 
			},{
				2*p1*u + 2*(3*Math.pow(rSq,2d)*k3*u + 2*(rSq)*k2*u + k1*u)*v + 2*p2*v,
				Math.pow(rSq,3d)*k3 + Math.pow(rSq,2d)*k2 + (rSq)*k1 + 2*p2*u + 2*(3*Math.pow(rSq,2d)*k3*v + 2*(rSq)*k2*v + k1*v)*v + 6*p1*v + 1
			}
		};
	}
	
	/** 
		Taken from the ImageJ source (http://imagej.nih.gov/ij/)
		modified
	*/
	public static byte[] getBilinearInterpolatedPixel(double x0, double y0, byte[][][] data) {
		int u0 = (int) Math.floor(x0);	//use floor to handle negative coordinates too
		int v0 = (int) Math.floor(y0);
		int width = data[0].length;
		int height = data.length;
		byte[] interp = new byte[3];
		
		if (u0>=0 && u0<width-1 && v0>= 0 && v0<height-1){
			double x = (x0-(double)u0);
			double y = (y0-(double)v0);
			for (int i =0;i<3;++i){
				interp[i] =(byte) (
						((double) data[v0][u0][i])*(1-x)*(1-y) 	/*f(0,0)(1-x)(1-y)*/
						+((double) data[v0+1][u0][i])*(1-y)*x	/*f(1,0)x(1-y)*/
						+((double) data[v0][u0+1][i])*(1-x)*y	/*f(0,1)(1-x)y*/
						+((double) data[v0+1][u0+1][i])*x*y	/*f(1,1)xy*/
						);
			}
		}
		return interp; /*Return zero for points outside the interpolable area*/
	}
	

	
		/*Taken from http://stackoverflow.com/questions/6524196/java-get-pixel-array-from-image
			Modified
		*/
   private static byte[][][] getPixelArray(BufferedImage image) {

		final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
		final int width = image.getWidth();
		final int height = image.getHeight();
		final boolean hasAlphaChannel = image.getAlphaRaster() != null;
		int pixelLength = 4;
		if (image.getAlphaRaster()!= null){
			pixelLength = 4;
		}else{
			pixelLength = 3;
		}
		byte[][][] result = new byte[height][width][3];
		for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
			result[row][col][0] = pixels[pixel + pixelLength-3]; //b
			result[row][col][1] = pixels[pixel + pixelLength-2]; //g
			result[row][col][2] = pixels[pixel + pixelLength-1]; //r
			col++;
			if (col == width) {
			   col = 0;
			   row++;
			}
		}
        return result;
	}
	
	/*Taken from http://stackoverflow.com/questions/6524196/java-get-pixel-array-from-image*/
   private static int[][] getPixelArrayInt(BufferedImage image) {

      final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
      final int width = image.getWidth();
      final int height = image.getHeight();
      final boolean hasAlphaChannel = image.getAlphaRaster() != null;

      int[][] result = new int[height][width];
      if (hasAlphaChannel) {
         final int pixelLength = 4;
         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
            int argb = 0;
            argb += (((int) pixels[pixel] & 0xff) << 24); // alpha
            argb += ((int) pixels[pixel + 1] & 0xff); // blue
            argb += (((int) pixels[pixel + 2] & 0xff) << 8); // green
            argb += (((int) pixels[pixel + 3] & 0xff) << 16); // red
            result[row][col] = argb;
            col++;
            if (col == width) {
               col = 0;
               row++;
            }
         }
      } else {
         final int pixelLength = 3;
         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
            int argb = 0;
            argb += -16777216; // 255 alpha
            argb += ((int) pixels[pixel] & 0xff); // blue
            argb += (((int) pixels[pixel + 1] & 0xff) << 8); // green
            argb += (((int) pixels[pixel + 2] & 0xff) << 16); // red
            result[row][col] = argb;
            col++;
            if (col == width) {
               col = 0;
               row++;
            }
         }
      }

      return result;
   }
}
