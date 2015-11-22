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
import java.util.Vector;
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

		JComponent cp = new JPanel();
		cp.add(loadImage("octaveTest/GOPR0093.JPG"));
		cp.add(loadImage("octaveTest/GOPR0099.JPG"));
		cp.setOpaque(true); //content panes must be opaque	
		setContentPane(cp);
		pack();
	        setVisible(true);
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
