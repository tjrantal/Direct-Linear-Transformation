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
	double scaleFactor;
	public ImagePanel(BufferedImage biIn){
		//Scale the image
		Image scale =  biIn.getScaledInstance(-1,250,Image.SCALE_SMOOTH);
		BufferedImage tempBi = new BufferedImage(scale.getWidth(null),scale.getHeight(null),biIn.getType());
		tempBi.getGraphics().drawImage(scale, 0, 0 , null);
		bi = tempBi;
		scaleFactor = ((double) bi.getWidth())/((double) biIn.getWidth());
		setPreferredSize(new Dimension(bi.getWidth(),bi.getHeight()));
		//paintImageToDraw();
	}
	
	public void paintImageToDraw(){
		repaint();
	}

	public void plotCoordinates(double[][] coordinates){
		plotCoordinates(coordinates,new Color(255,0,0));
	}
	public void plotCoordinates(double[][] coordinates,Color color){
		Graphics2D g2 = bi.createGraphics();
		g2.setColor(color);
		g2.setStroke(new BasicStroke());
		int x,y;
		for (int i = 0; i<coordinates.length;++i){
			x = (int) (coordinates[i][0]*scaleFactor);
			y = (int) (coordinates[i][1]*scaleFactor);
			//System.out.println("Plotting "+x+"\t"+y);
			g2.drawOval(x-2,y-2,4,4);
		}
		g2.dispose();
		paintImageToDraw();
	}

	public void paint(Graphics g) {
		g.drawImage(bi,0,0,null);
	}
}
