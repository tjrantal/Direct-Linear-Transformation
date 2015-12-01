package timo.home.undristort;

import Jama.*;
import java.awt.image.*;

public class Undistort{
	BufferedImage ubi;	//Undistorted buffered image
	/*Constructor*/
	public Undistort(BufferedImage bi, double scaleFactor, Matrix K, Matrix d){
		byte[][][] pixArr = getPixelArray(bi);
	
	}
	
	/** This method is from Chapter 16 of "Digital Image Processing:
	An Algorithmic Introduction Using Java" by Burger and Burge
	(http://www.imagingbook.com/). 
		Taken from the ImageJ source (http://imagej.nih.gov/ij/)
	*/
	public static double getBicubicInterpolatedPixel(double x0, double y0, double[][] data) {
		int u0 = (int) Math.floor(x0);	//use floor to handle negative coordinates too
		int v0 = (int) Math.floor(y0);
		int width = data.length;
		int height = data[0].length;
		if (u0<1 || u0>width-3 || v0< 1 || v0>height-3){
			if ((u0 == 0 || u0 < width-1) && (v0 == 0 || v0 < height-1)){ /*Use bilinear interpolation http://en.wikipedia.org/wiki/Bilinear_interpolation*/
				double x = (x0-(double)u0);
				double y = (y0-(double)v0);
				return data[u0][v0]*(1-x)*(1-y) 	/*f(0,0)(1-x)(1-y)*/
						+data[u0+1][v0]*(1-y)*x	/*f(1,0)x(1-y)*/
						+data[u0][v0+1]*(1-x)*y	/*f(0,1)(1-x)y*/
						+data[u0+1][v0+1]*x*y;	/*f(1,1)xy*/
			}
			return 0; /*Return zero for points outside the interpolable area*/
		}
		double q = 0;
		for (int j = 0; j < 4; ++j) {
			int v = v0 - 1 + j;
			double p = 0;
			for (int i = 0; i < 4; ++i) {
				int u = u0 - 1 + i;
				p = p + data[u][v] * cubic(x0 - u);
			}
			q = q + p * cubic(y0 - v);
		}
		return q;
	}
	
	public static final double cubic(double x) {
		final double a = 0.5; // Catmull-Rom interpolation
		if (x < 0.0) x = -x;
		double z = 0.0;
		if (x < 1.0) 
			z = x*x*(x*(-a+2.0) + (a-3.0)) + 1.0;
		else if (x < 2.0) 
			z = -a*x*x*x + 5.0*a*x*x - 8.0*a*x + 4.0*a;
		return z;
	}
	
		/*Taken from http://stackoverflow.com/questions/6524196/java-get-pixel-array-from-image
			Modified
		*/
   private static byte[][][] getPixelArray(BufferedImage image) {

      final byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();
      final int width = image.getWidth();
      final int height = image.getHeight();
      final boolean hasAlphaChannel = image.getAlphaRaster() != null;

      byte[][][] result = new int[height][width][3];
      if (hasAlphaChannel) {
         final int pixelLength = 4;
         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
            results[row][col][0] = pixels[pixel + 1]; //b
            results[row][col][1] = pixels[pixel + 2]; //g
            results[row][col][2] = pixels[pixel + 3]; //r
            col++;
            if (col == width) {
               col = 0;
               row++;
            }
         }
      } else {
         final int pixelLength = 3;
         for (int pixel = 0, row = 0, col = 0; pixel < pixels.length; pixel += pixelLength) {
            results[row][col][0] = pixels[pixel]; //b
            results[row][col][1] = pixels[pixel + 1]; //g
            results[row][col][2] = pixels[pixel + 2]; //r
            col++;
            if (col == width) {
               col = 0;
               row++;
            }
         }
      }
      return result;
   }
	
	/*Taken from http://stackoverflow.com/questions/6524196/java-get-pixel-array-from-image*/
   private static int[][] getPixelArray(BufferedImage image) {

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
   
   /**
   *Taken from http://stackoverflow.com/questions/5617016/how-do-i-copy-a-2-dimensional-array-in-java
	 * Clones the provided array
	 * 
	 * @param src
	 * @return a new clone of the provided array
	 */
	public static int[][] cloneArray(int[][] src) {
		 int length = src.length;
		 int[][] target = new int[length][src[0].length];
		 for (int i = 0; i < length; i++) {
		     System.arraycopy(src[i], 0, target[i], 0, src[i].length);
		 }
		 return target;
	}
}
