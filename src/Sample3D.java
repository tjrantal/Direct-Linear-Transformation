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

import Jama.*;
import dlt.*;
import java.util.Vector;

public class Sample3D{
	 public static void main(String[] args) {
	 			//Laitetaan koordinaatit paikoilleen
		double[][] global = {{0,0,0},{10,0,0},{0,10,0},{10,10,0},{0,0,10},{10,0,10}}; //Global coordinates of the calibration object
		double[][][] calib={
											{{823,591},{833,351},{967,585},{981,330},{616,573},{628,322}},
											{{1300,647},{1200,484},{1429,492},{1320,327},{1143,654},{1046,480}}
										};	//Digitized calibration object points
		double[][] digitizedUnknownPoint = {{855,484},{1292,470}};	//Digitized unknown point N,2 Matrix
		System.out.println("DLT3D");
		DLT3D dlt3d = new DLT3D(global,calib);
		Matrix coordinates = dlt3d.scaleCoordinates(digitizedUnknownPoint);
		
		/*Print out results*/
		Vector<Matrix> coeffs = dlt3d.getCurrentDltCoefficients(); 
		String resultString = "Coefficients\n";
		for (int c = 0; c< coeffs.size();++c){
			resultString+="\tCam"+c+"\n";
			for (int i = 0; i< coeffs.get(c).getRowDimension();++i){
				resultString+="\t\t"+i+": "+coeffs.get(c).get(i,0)+"\n"; 
			}
		}
		System.out.println(resultString);
		System.out.println("X "+coordinates.get(0,0)+" Y "+coordinates.get(1,0)+" "+coordinates.get(2,0));		
	}
}
