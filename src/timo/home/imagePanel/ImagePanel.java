package timo.home.imagePanel;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.Vector;
import java.awt.BasicStroke;
import java.awt.image.*;
import java.awt.image.WritableRaster;
import java.text.DecimalFormat;	//For rounding text

public class ImagePanel extends JPanel{
	BufferedImage bi;
	public ImagePanel(BufferedImage biIn){
		//Scale the image
		Image scale =  biIn.getScaledInstance(-1,250,Image.SCALE_SMOOTH);
		BufferedImage tempBi = new BufferedImage(scale.getWidth(null),scale.getHeight(null),biIn.getType());
		tempBi.getGraphics().drawImage(scale, 0, 0 , null);
		bi = tempBi;
		setPreferredSize(new Dimension(bi.getWidth(),bi.getHeight()));
		paintImageToDraw();
	}
	
	public void paintImageToDraw(){
		repaint();
	}

	public void paint(Graphics g) {
		g.drawImage(bi,0,0,null);
	}
}
