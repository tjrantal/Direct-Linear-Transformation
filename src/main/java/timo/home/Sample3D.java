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
	Direct linear transformation (DLT) sample, written by Timo Rantalainen tjrantal@gmail.com.
	Depends on Jama library http://math.nist.gov/javanumerics/jama/
	compile javac Sample3D.java
	execute java Sample3D
*/

package timo.home;

import Jama.*;
import timo.home.dlt.*;
import timo.home.imagePanel.*;
import timo.home.utils.*;
import timo.home.calibration.*;
import timo.home.refine.*;
import timo.home.undistort.*;

import java.util.*;
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
import javax.swing.*;

public class Sample3D extends JFrame{
	/*Constructor*/
	public Sample3D(String title){
		super(title);
		setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		//Read digitized data
		CSVReader[] cr = new CSVReader[3];
	 	cr[0] = new CSVReader("octaveTest/sampleFigs/GOPR0093.txt","\t");
		cr[1] = new CSVReader("octaveTest/sampleFigs/GOPR0099.txt","\t");
		cr[2] = new CSVReader("octaveTest/calibrationObject.xls","\t");
		
		//Add images
		JComponent cp = new JPanel();
		ImagePanel[] ip = new ImagePanel[4];
		//Scale images down to 1/5th. My laptop couldn't handle the full image...
		ImagePanel tempI = loadImage("octaveTest/sampleFigs/GOPR0093.JPG");
		ip[0] = new ImagePanel(tempI.scaleImage(tempI.getOrig(),0.2));
		ip[0].setOpaque(true);
		tempI = loadImage("octaveTest/sampleFigs/GOPR0099.JPG");
		ip[1] = new ImagePanel(tempI.scaleImage(tempI.getOrig(),0.2));
		ip[1].setOpaque(true);
		
		double[][] calibrationObject = getCoords(cr[2]);
		/*
		System.out.println("");
		for (int i = 0; i<calibrationObject.length;++i){
			System.out.print(i+"\t");
			for (int j = 0; j<calibrationObject[i].length;++j){
				System.out.print(calibrationObject[i][j]+"\t");
			}
			System.out.println("");
			
		}
		*/
		LensCalibration lc = new LensCalibration(calibrationObject,ip[0].origWidth/2,ip[0].origHeight/2);
		double[][][] undistortedCoord = new double[2][][];
		double[][][] undistortedUnknown = new double[2][][];
		Undistort[] ud = new Undistort[2];
		Matrix[][] Rt = new Matrix[2][2];
		double[][] apacheCoeffs = new double[2][];
		for (int i = 0;i<2;++i){
			
			
			double[][] digitized = getCoords(cr[i]);
			double[][] confirm = new double[digitized.length][2];
			System.out.println(String.format("Cam %02d length digitized %d",i,digitized.length));
			//Scale the coordinates down by 5. My laptop couldn't handle the full image...
			for (int r = 0;r<digitized.length;++r){
				//System.out.print(r+"\t");
				for (int c = 4;c<digitized[r].length;++c){
					//System.out.print(digitized[r][c]+"\t");
					digitized[r][c-4] = digitized[r][c]/5d;	//The .txt files contain multiple rows -> pop the X and Y into the correct columns and scale down by five..
				}
				
				//Set 50% of coordinates to null
				if (r % 2 > 0){
					confirm[r][0] = digitized[r][0];
					confirm[r][1] = digitized[r][1]; 
					digitized[r] = null;
				}else{
					confirm[r] = null;
				}
				//System.out.println("");	
			}
			
			ArrayList<Matrix> KRt = lc.calibrate(digitized);
			
			if(false){
				System.out.println("K");
				KRt.get(0).print(3,3);
				System.out.println("R");
				KRt.get(1).print(3,3);
				System.out.println("t");
				KRt.get(2).print(3,3);
				System.out.println("");
			}
			
			//Refine the calibration
			RefineParameters rp = new RefineParameters(calibrationObject,digitized,KRt.get(0),new double[]{0d,0d,0d,0d,0d},KRt.get(1),KRt.get(2));
			ArrayList<Matrix> KdRt = rp.getCalibration();
			Rt[i][0]=KdRt.get(2);
			Rt[i][1]=KdRt.get(3);
			
			if(false){
				System.out.println("refined K");
				KdRt.get(0).print(3,3);
				System.out.println("d");
				KdRt.get(1).print(3,3);
				System.out.println("R");
				KdRt.get(2).print(3,3);
				System.out.println("t");
				KdRt.get(3).print(3,3);
			}
			
			
			ud[i] = new Undistort(ip[i], KdRt.get(0), KdRt.get(1));
			ip[i+2] = new ImagePanel(ud[i].ubi);
			ip[i].plotCoordinates(digitized);
			undistortedCoord[i] = ud[i].undistortCoordinates(digitized);
			ip[i+2].plotCoordinates(undistortedCoord[i]);
			undistortedUnknown[i] = ud[i].undistortCoordinates(confirm);
			apacheCoeffs[i] = DLT3D_apache_math3.getDltCoefficients(calibrationObject,undistortedCoord[i]);	//Compute Apache DLT


			System.out.println("t");
			KRt.get(2).print(3,3);
			KdRt.get(3).print(3,3);			
			System.out.println("");
			double[] camPos = DLT3D_apache_math3.getCameraPosition(apacheCoeffs[i]);
			for (double a:camPos){
				System.out.println(String.format("%.03f",a));
			}


			double[][] projectedPoints =  ud[i].projectKnownPoints(calibrationObject, KdRt.get(2), KdRt.get(3));
			double[][] apacheProjectedPoints = new double[calibrationObject.length][];
			int cnt = 0;
			System.out.println("Projected object "+i);
			for (double[] point:calibrationObject){
				apacheProjectedPoints[cnt] = DLT3D_apache_math3.project3D(apacheCoeffs[i], point);
				System.out.println(String.format("p%02d\t%.2f\t%.2f\t%.2f\t%.2f",cnt,projectedPoints[cnt][0],projectedPoints[cnt][1],apacheProjectedPoints[cnt][0],apacheProjectedPoints[cnt][1]));
				++cnt;
			}

			ip[i+2].plotCoordinates(projectedPoints, new Color(0,255,0));
			ip[i+2].plotCoordinates(apacheProjectedPoints, new Color(255,0,255));
			
			//ip[i+2].paintImageToDraw();
			//ip[i+2].setOpaque(true);
			
			cp.add(ip[i]);
			cp.add(ip[i+2]);
			
			//DLT based on undistorted coordinates
			
		}
		DLT3D dlt3d = new DLT3D(calibrationObject,undistortedCoord);
		Vector<Matrix> dltCoeffs= dlt3d.getCurrentDltCoefficients();
		
		float[] dltCoeffsToJOCL = new float[dltCoeffs.size()*11];	//number_of_cams x 11 coefficients
		for (int c = 0; c<dltCoeffs.size();++c){
			System.out.println(String.format("Cam %d\tJama\tGradle",c));
			
			for (int i = 0; i<apacheCoeffs[c].length;++i){
				System.out.println(String.format("coeff%02d\t%.2f\t%.2f",c,dltCoeffs.get(c).get(i,0),apacheCoeffs[c][i]));
				dltCoeffsToJOCL[c*11+i] = (float) apacheCoeffs[c][i];
			}
		}

		//Get 3D coords, and project to undistorted coords. UnknownPoints need to be undistorted as well
		double[][] unknownUndistorted3D = new double[undistortedUnknown[0].length][3];
		Matrix tempC;
		float[] camCoordinatesToJOCL = new float[dltCoeffs.size()*undistortedUnknown[0].length*2]; //number_of_digitised_points x camNo x 2 coordinates
		for (int i = 0;i<unknownUndistorted3D.length;++i){	//Digitised points
			for (int c = 0;c<dltCoeffs.size();++c){//Camera views
				if (undistortedUnknown[c][i] != null){
					camCoordinatesToJOCL[i*2*dltCoeffs.size()+2*c+0] = (float) undistortedUnknown[c][i][0];
					camCoordinatesToJOCL[i*2*dltCoeffs.size()+2*c+1] = (float) undistortedUnknown[c][i][1];
				}else{
					camCoordinatesToJOCL[i*2*dltCoeffs.size()+2*c+0] = 0f;
					camCoordinatesToJOCL[i*2*dltCoeffs.size()+2*c+1] = 0f;
				}
			}
		}
		//Run 2D to 3D with JOCL
		float[] jocl3d = JOCL_views_to_3D.testJOCL(2,undistortedUnknown[0].length,camCoordinatesToJOCL,dltCoeffsToJOCL);
		
		for (int i =0;i<unknownUndistorted3D.length;++i){
			if (undistortedUnknown[0][i] != null){
				tempC = dlt3d.scaleCoordinates(new double[][]{undistortedUnknown[0][i],undistortedUnknown[1][i]});
				unknownUndistorted3D[i] = new double[]{tempC.get(0,0),tempC.get(1,0),tempC.get(2,0)};
				double[] apache3D = DLT3D_apache_math3.scaleCoordinates(apacheCoeffs,new double[][]{undistortedUnknown[0][i],undistortedUnknown[1][i]});
				System.out.println(String.format("Jama %d x %.1f y %.1f z %.1f",i,tempC.get(0,0),tempC.get(1,0),tempC.get(2,0)));
				System.out.println(String.format("Apache %d x %.1f y %.1f z %.1f",i,apache3D[0],apache3D[1],apache3D[2]));
			}else{
				unknownUndistorted3D[i] = null;
			}
		}	
		
		ip[2].plotCoordinates(ud[0].projectKnownPoints(unknownUndistorted3D, Rt[0][0], Rt[0][1]), new Color(0,0,255));
		ip[3].plotCoordinates(ud[1].projectKnownPoints(unknownUndistorted3D, Rt[1][0], Rt[1][1]), new Color(0,0,255));
		//Matrix coordinates = dlt3d.scaleCoordinates(digitizedUnknownPoint);
		cp.setOpaque(true); // must be opaque	
		setContentPane(cp);
		pack();
		setSize(new Dimension(680,550));
	   setVisible(true);

	   
	}

	private double[][] getCoords(CSVReader cr){
		double[][] temp = new double[cr.data.size()-1][cr.data.get(0).size()-1];
		for (int i = 1; i< cr.data.size();++i){
			for (int j = 1; j< cr.data.get(i).size();++j){
				temp[i-1][j-1] = Double.parseDouble(cr.data.get(i).get(j));
			}
		}
		return temp;
	} 
	private ImagePanel loadImage(String fileName){
		BufferedImage img = null;
		try {
		    img = ImageIO.read(new File(fileName));
		} catch (IOException e) {System.out.println("Couldn't read image");}
		ImagePanel photo = new ImagePanel(img);
		photo.setOpaque(true);
		return photo;
	}
	/*Main to get things going*/	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable(){
      			public void run(){
        			Sample3D mainClass = new Sample3D("Sample 3d");
        		}
    		});
  	}
	
}
