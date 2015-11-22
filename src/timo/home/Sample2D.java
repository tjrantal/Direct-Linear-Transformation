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
	compile javac Sample2D.java
	execute java Sample2D
*/

package timo.home;

import Jama.*;
import timo.home.dlt.*;

public class Sample2D{
	 public static void main(String[] args) {
		double[][] global =  {{0,0}, {29.6, 0}, {0, 21}, {29.6, 21}};  //Global coordinates of the calibration object
		double[][] calib = {{900,724},{2320, 575},{955, 1711},{ 2610, 1501}};  //Digitized calibration object points
		double[]  digitizedUnknownPoint = {1245.0,919.0};	//Digitized unknown point
		System.out.println("DLT2D");
		DLT2D dlt2d = new DLT2D(global,calib);
		Matrix coordinates = dlt2d.scaleCoordinates(digitizedUnknownPoint);
		
		Matrix coeffs = dlt2d.getCurrentDltCoefficients();
		String resultString = "Coefficients\n";
		for (int i = 0; i< coeffs.getRowDimension();++i){
			resultString+="\t"+i+": "+coeffs.get(i,0)+"\n"; 
		}
		System.out.println(resultString);
		System.out.println("X "+coordinates.get(0,0)+" Y "+coordinates.get(1,0));
	}
}
