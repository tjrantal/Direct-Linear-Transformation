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
	
	public double[][] undistortCoordinates(double[][] coordinates){
		double[][] undistorted = new double[coordinates.length][coordinates[0].length];
		double x_u,y_u,radius,radial;
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
		for (int i = 0;i<coordinates.length;++i){
			x_u = xx_up[0][i];
			y_u = xx_up[1][i];
			radius = Math.sqrt(x_u*x_u + y_u*y_u);
			radial = 1d/(1d + k1*radius*radius + k2*radius*radius*radius*radius  + k3*radius*radius*radius*radius*radius*radius);
			xxdTemp[0][i] = radial*x_u - 2*p1*x_u*y_u - p2*(radius*radius + 2*x_u*x_u);
			xxdTemp[1][i] = radial*y_u - p1*(radius*radius + 2*y_u*y_u) - 2*p2*x_u*y_u;
			xxdTemp[2][i] = 1d;
		}
		Matrix xx_dim = K.times(new Matrix(xxdTemp));
		double[][] xx_dimp = xx_dim.getArray();
		for (int i = 0;i<coordinates.length;++i){
				undistorted[i][0] = xx_dimp[0][i];
				undistorted[i][1] = xx_dimp[1][i];
		}
		return undistorted;
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
