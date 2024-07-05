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
	BufferedImage orig;
	double scaleFactor;
	public int origWidth;
	public int origHeight;
	public ImagePanel(BufferedImage biIn){
		//Scale the image
		orig = biIn;
		origWidth = biIn.getWidth();
		origHeight = biIn.getHeight();
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
			if (coordinates[i] != null){
				x = (int) (coordinates[i][0]*scaleFactor);
				y = (int) (coordinates[i][1]*scaleFactor);
				System.out.println("Plotting "+x+"\t"+y);
				g2.drawOval(x-1,y-1,3,3);
			}
		}
		g2.dispose();
		paintImageToDraw();
	}

	public void paint(Graphics g) {
		g.drawImage(bi,0,0,null);
	}
	
	public BufferedImage getBI(){
		return bi;
	}
	
	public BufferedImage getOrig(){
		return orig;
	}
	
	public void setBI(BufferedImage bi){
		this.bi = bi;
	}
	
	public double getScaleFactor(){
		return scaleFactor;
	}
	
	public BufferedImage scaleImage(BufferedImage biIn,double relativeSize){
		Image scale =  biIn.getScaledInstance((int) (biIn.getWidth()*relativeSize),-1,Image.SCALE_SMOOTH);
		BufferedImage tempBi = new BufferedImage(scale.getWidth(null),scale.getHeight(null),biIn.getType());
		tempBi.getGraphics().drawImage(scale, 0, 0 , null);
		return tempBi;
	}
	
	public BufferedImage createImageFromBytes(byte[][][] dataSource){
		int w = dataSource[0].length;
		int h = dataSource.length;
		int[] pix = new int[w*h];
		for (int r =0;r<h;++r){
			for (int c =0;c<w;++c){
				pix[c+r*w] = (255 << 24) | ((((int) dataSource[r][c][2]) & 0xff) << 16) | ((((int) dataSource[r][c][1]) & 0xff) << 8) | ((((int) dataSource[r][c][0] & 0xff)));
				//pix[c+r*w] = (255 << 24) | 255 << 16;
			}
		}
		Image img = createImage(new MemoryImageSource(w, h, pix, 0, w));
		BufferedImage temp = new BufferedImage(img.getWidth(null),img.getHeight(null),BufferedImage.TYPE_INT_ARGB);
		temp.getGraphics().drawImage(img, 0, 0 , null);
		return temp;
	}
}
