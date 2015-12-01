package timo.home.undristort;

import Jama.*;
import java.awt.image.*;

public class Undistort{
	BufferedImage ubi;	//Undistorted buffered image
	/*Constructor*/
	public Undistort(BufferedImage bi, double scaleFactor, Matrix K, Matrix d){
		byte[][][] pixArr = getPixelArray(bi);
	
	}

	/** 
		Taken from the ImageJ source (http://imagej.nih.gov/ij/)
		modified
	*/
	public static byte[] getBilinearInterpolatedPixel(double x0, double y0, byte[][][] data) {
		int u0 = (int) Math.floor(x0);	//use floor to handle negative coordinates too
		int v0 = (int) Math.floor(y0);
		int width = data.length;
		int height = data[0].length;
		byte[] interp = new byte[3];
		
		if (u0>=0 && u0<width-1 && v0>= 0 && v0<height-1){
			double x = (x0-(double)u0);
			double y = (y0-(double)v0);
			for (int i =0;i<3;++i){
				interp[i] =(byte) (
						((double) data[u0][v0][i])*(1-x)*(1-y) 	/*f(0,0)(1-x)(1-y)*/
						+((double) data[u0+1][v0][i])*(1-y)*x	/*f(1,0)x(1-y)*/
						+((double) data[u0][v0+1][i])*(1-x)*y	/*f(0,1)(1-x)y*/
						+((double) data[u0+1][v0+1][i])*x*y	/*f(1,1)xy*/
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
}
