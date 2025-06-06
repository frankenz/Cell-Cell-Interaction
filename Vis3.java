import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.imageio.*;
import java.awt.image.*;
import javax.swing.filechooser.*;
import java.util.*;
import java.io.*;

public class Vis3 extends JApplet implements ActionListener{
	public static String version="21st June 2012";
	final static Color bg = Color.white;
	final static Color fg = Color.black;
	Dimension totalSize;
	public BufferedImage img=null; 
	public int size=1000;
	int numParticles=0;
	int[][] OriginalMatrix;
	int[][] imgBuffer;
	int mag=2; // Magnification factor 
	public prosCA3 ca;
	int counter=0; // Timesteps counter
	int defaultView=0; // Cells view
	JComboBox views=null;
	JButton button;
	int timestep=0;
	boolean movie=false;
	public boolean play=false;
	JFileChooser fc;
	JLabel label;
	JCheckBox cb;
	String rootDir=".";


	public Vis3 (JComboBox v, JButton b, JLabel l, JCheckBox _cb)
	{
		fc = new JFileChooser();
		button = b;
		label = l;
		views=v;
		cb=_cb;
		ca = new prosCA3();
		imgBuffer = new int [ca.size][ca.size];
		size=ca.size;
		mag = ca.mag;
		img = new BufferedImage (4*ca.size*ca.mag,2*ca.size*ca.mag,BufferedImage.TYPE_INT_RGB);
	}
	
	public final BufferedImage scale(double scale, BufferedImage srcImg) 
	{ 
		if (scale == 1)  return srcImg; 
   		AffineTransformOp op = new AffineTransformOp( AffineTransform.getScaleInstance( scale, scale), null); 
        	return op.filter(srcImg, null); 
	}

	public void paint(Graphics g) {
        	Graphics2D g2 = (Graphics2D) g;
			BufferedImage img2=scale (mag,img);
			g2.drawImage(img2,null,null);
			if (ca.size!=size) {
				imgBuffer = new int [ca.size][ca.size];
				img = new BufferedImage (3*ca.size*ca.mag,2*ca.size*ca.mag,BufferedImage.TYPE_INT_RGB);
				ca.size=size;
			}       
	}



    public void nextTimeStep ()
	{
		ca.nextTimeStep();
		if (defaultView==0) paintCells();
		timestep++;
	}
    

	public void paintCells()
	{
		int [][] lattice = ca.getCells();
        
		
		BufferedImage all = new BufferedImage (ca.size*3,ca.size*2,BufferedImage.TYPE_INT_RGB);  //NEW FROM HERE//////////////////////////////////////////
                
		BufferedImage c= getCells();
		BufferedImage b=getBM();
		BufferedImage t=getTGFB();
		BufferedImage m=getMDE();

 		for (int i=0;i<ca.size;i++)
			for (int j=0;j<size;j++) { 
                               
            /*    if ((lattice[i][j]==3) && (ca.MDEProduction[i][j]<0.33)) img.setRGB(i+size+size,j+size,Color.green.getRGB());
                else if ((lattice[i][j]==3) && (ca.MDEProduction[i][j]>=0.33)&& (ca.MDEProduction[i][j]<0.66)) img.setRGB(i+size+size,j+size,Color.yellow.getRGB());
                else if ((lattice[i][j]==3) && (ca.MDEProduction[i][j]>=0.66)) img.setRGB(i+size+size,j+size,Color.red.getRGB());
				else if (lattice[i][j]==3)img.setRGB(i+size+size,j+size,Color.white.getRGB());
                else img.setRGB(i+size+size,j+size,c.getRGB(i,j));   */
                
                img.setRGB(i+size+size,j+size,Color.white.getRGB());
                if (lattice[i][j]==1) img.setRGB(i+size+size,j+size,Color.black.getRGB()); // NEW Basel cells
				if (lattice[i][j]==2) img.setRGB(i+size+size,j+size,Color.gray.getRGB()); // NEW Lumen cells
                
                if (lattice[i][j]==3){
				if (ca.MDEProduction[i][j]<=0.01f) img.setRGB(i+size+size,j+size,Color.white.getRGB());
				else if (ca.MDEProduction[i][j]<=0.15) {
					Color darkBlue = new Color(0, 0, 128);  
                    img.setRGB(i+size+size,j+size,darkBlue.getRGB());
                }
                else if (ca.MDEProduction[i][j]<=0.2) img.setRGB(i+size+size,j+size,Color.blue.getRGB());
				else if (ca.MDEProduction[i][j]<=0.25) img.setRGB(i+size+size,j+size,Color.cyan.getRGB());
				else if (ca.MDEProduction[i][j]<=0.3) {
					Color darkCyan = new Color(0, 128, 128);  
                    img.setRGB(i+size+size,j+size,darkCyan.getRGB());
                }
				else if (ca.MDEProduction[i][j]<=0.5) img.setRGB(i+size+size,j+size,Color.green.getRGB());
				else if (ca.MDEProduction[i][j]<=0.6) {
					Color limeGreen = new Color(170, 255, 0);  
                    img.setRGB(i+size+size,j+size,limeGreen.getRGB());
                }
				else if (ca.MDEProduction[i][j]<=0.65) img.setRGB(i+size+size,j+size,Color.yellow.getRGB());
				else if (ca.MDEProduction[i][j]<=0.7) img.setRGB(i+size+size,j+size,Color.orange.getRGB());
				else if (ca.MDEProduction[i][j]<=1.0) img.setRGB(i+size+size,j+size,Color.red.getRGB()); 
                else img.setRGB(i+size+size,j+size,c.getRGB(i,j));
                                        
                }
              /*  if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) img.setRGB(i+size+size,j+size,Color.white.getRGB()); // border  */

                
            /*    if ((lattice[i][j]==3) && (ca.TGFBProduction[i][j]<0.33)) img.setRGB(i+size+size,j,Color.green.getRGB());
                else if ((lattice[i][j]==3) && (ca.TGFBProduction[i][j]>=0.33)&& (ca.TGFBProduction[i][j]<0.66)) img.setRGB(i+size+size,j,Color.yellow.getRGB());
                else if ((lattice[i][j]==3) && (ca.TGFBProduction[i][j]>=0.66)) img.setRGB(i+size+size,j,Color.red.getRGB());
				else if (lattice[i][j]==3) img.setRGB(i+size+size,j,Color.white.getRGB());
                else img.setRGB (i+size+size,j,c.getRGB(i,j));
             */

                img.setRGB(i+size+size,j,Color.white.getRGB());
                if (lattice[i][j]==1) img.setRGB(i+size+size,j,Color.black.getRGB()); // NEW Basel cells
				if (lattice[i][j]==2) img.setRGB(i+size+size,j,Color.gray.getRGB()); // NEW Lumen cells

                if (lattice[i][j]==3){
                    if (ca.TGFBProduction[i][j]<=0.01f) img.setRGB(i+size+size,j,Color.white.getRGB());
                    else if (ca.TGFBProduction[i][j]<=0.15) {
                        Color darkBlue = new Color(0, 0, 128);  
                        img.setRGB(i+size+size,j,darkBlue.getRGB());
                    }
                    else if (ca.TGFBProduction[i][j]<=0.2) img.setRGB(i+size+size,j,Color.blue.getRGB());
                    else if (ca.TGFBProduction[i][j]<=0.25) img.setRGB(i+size+size,j,Color.cyan.getRGB());
                    else if (ca.TGFBProduction[i][j]<=0.3) {
                        Color darkCyan = new Color(0, 128, 128);  
                        img.setRGB(i+size+size,j,darkCyan.getRGB());
                    }
                    else if (ca.TGFBProduction[i][j]<=0.5) img.setRGB(i+size+size,j,Color.green.getRGB());
                    else if (ca.TGFBProduction[i][j]<=0.6) {
                        Color limeGreen = new Color(170, 255, 0);  
                        img.setRGB(i+size+size,j,limeGreen.getRGB());
                    }
                    else if (ca.TGFBProduction[i][j]<=0.65) img.setRGB(i+size+size,j,Color.yellow.getRGB());
                    else if (ca.TGFBProduction[i][j]<=0.7) img.setRGB(i+size+size,j,Color.orange.getRGB());
                    else if (ca.TGFBProduction[i][j]<=1.0) img.setRGB(i+size+size,j,Color.red.getRGB());
                    else img.setRGB (i+size+size,j,c.getRGB(i,j));
                    
                    
                }                    
             /*   if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) img.setRGB(i+size+size,j,Color.white.getRGB()); // border */

img.setRGB(i+size+size+size,j,Color.white.getRGB());
                if (lattice[i][j]==1) img.setRGB(i+size+size+size,j,Color.black.getRGB()); // NEW Basel cells
				if (lattice[i][j]==2) img.setRGB(i+size+size+size,j,Color.gray.getRGB()); // NEW Lumen cells
				if ( (lattice[i][j]==6) && (ca.stromalCellTGFBproductionRate > 0.19) ) {
                    Color darkRed = new Color(170, 0, 0);  
                    img.setRGB(i+size+size+size,j,darkRed.getRGB());
                } // Act Stromal
                if ( (lattice[i][j]==6) && (ca.stromalCellTGFBproductionRate > 0.1) && (ca.stromalCellTGFBproductionRate < 0.13) ) {
                    Color darkGreen = new Color(0, 128, 0);  
                    img.setRGB(i+size+size+size,j,darkGreen.getRGB());
                } // Act Stromal
				if ( (lattice[i][j]==6) && (ca.stromalCellTGFBproductionRate < 0.01) ) img.setRGB(i+size+size+size,j,Color.blue.getRGB()); // Act Stromal



                all.setRGB(i,j,c.getRGB(i,j));
				all.setRGB(i+size,j,t.getRGB(i,j));
				all.setRGB (i,j+size,b.getRGB(i,j));
				all.setRGB (i+size,j+size,m.getRGB(i,j));
                img.setRGB(i+size,j,t.getRGB(i,j));
				img.setRGB (i,j+size,b.getRGB(i,j));
				img.setRGB (i+size,j+size,m.getRGB(i,j));
                img.setRGB(i,j,c.getRGB(i,j));
                
                // GET DUCTS TO ALL TO t and m 
                
                if (lattice[i][j]==1) img.setRGB(i+size,j+size,Color.black.getRGB()); // NEW Basel cells
				if (lattice[i][j]==2) img.setRGB(i+size,j+size,Color.gray.getRGB()); // NEW Lumen cells                
                
                if (lattice[i][j]==1) img.setRGB(i+size,j,Color.black.getRGB()); // NEW Basel cells
				if (lattice[i][j]==2) img.setRGB(i+size,j,Color.gray.getRGB()); // NEW Lumen cells                
                
    
                }
        
/*		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				int val=0;
				if (lattice[i][j]==0) val=Color.white.getRGB();
				else if (lattice[i][j]==1) val=Color.blue.getRGB(); // Basel cells
				else if (lattice[i][j]==2) val=Color.green.getRGB(); // Lumen cells
				else if (lattice[i][j]==3) val=Color.red.getRGB(); // Cancer (LU ab)
				else if (lattice[i][j]==4) {
					Color brown = new Color(255, 128, 0); 
					val=brown.getRGB();// Stromal
				}
				else if (lattice[i][j]==5) val=Color.black.getRGB();// Monocytes	
				img.setRGB(i,j,val);
			} */
    
        
                                        //NEW TO HERE ////////////////////////////////////////////////////////////////////////////////////////////
        
		repaint();
		if (movie) try {
            
            if ((counter%1==0) && (ca.totalBirths > 0)){   //NEW
                
			File dir = new File (rootDir+"/"+ca.OutName+"/"+ca.OutName+"Images"); // NEW
			dir.mkdir ();
			File fileAll = new File (rootDir+"/"+ca.OutName+"/"+ca.OutName+"Images/"+ca.OutName+"T"+counter+".png");   // NEW
			ImageIO.write(scale(mag,img),"png",fileAll);
            } //NEW
		} catch (Exception e) {
			e.printStackTrace();
		}	

		counter++;
	}

    
    
    public BufferedImage getCells()
	{
		int [][] lattice = ca.getCells();
		BufferedImage result = new BufferedImage (ca.size,ca.size,BufferedImage.TYPE_INT_RGB);
		
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				int val=0;
				if (lattice[i][j]==0) val=Color.white.getRGB();
				else if (lattice[i][j]==1) val=Color.black.getRGB(); // NEW Basel cells
				else if (lattice[i][j]==2) val=Color.gray.getRGB(); // NEW Lumen cells
                else if (lattice[i][j]==3) {
					Color brown = new Color(102, 0, 0); 
					val=brown.getRGB();
				} // Cancer (LU ab)
				else if (lattice[i][j]==4) {
					Color magenta = new Color(170, 0, 255); 
					val=magenta.getRGB();
				}// NEW Stromal
				else if (lattice[i][j]==5) {
					Color pink = new Color(255, 0, 128); 
					val=pink.getRGB();
				}// Monocytes
                		else if (lattice[i][j]==6) val=Color.yellow.getRGB(); // Active Stroma as now is stromalCellTGFBproductionRate
                /*  if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) val=Color.white.getRGB(); // border  */        
                result.setRGB(i,j,val);  
			}
		return result;
	}

	public BufferedImage getTGFB ()
	{
		BufferedImage result = new BufferedImage (ca.size,ca.size,BufferedImage.TYPE_INT_RGB);
        float[][] lattice = ca.getTGFB();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.15) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.25) c=Color.cyan;
                else if (lattice[i][j]<=0.3) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.8) c=Color.yellow;
				else if (lattice[i][j]<=0.995) c=Color.orange;
                else if (lattice[i][j]<=0.998) c= new Color(255, 69, 0); // redOrange
 				else if (lattice[i][j]<=1.0) c=Color.red;
                else System.out.println ("Error here");
              /*  if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                result.setRGB(i,j,c.getRGB());               
			}
		return result;
	}

	public void paintTGFB ()
	{
		float[][] lattice = ca.getTGFB();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.15) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.25) c=Color.cyan;
                else if (lattice[i][j]<=0.3) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.8) c=Color.yellow;
				else if (lattice[i][j]<=0.995) c=Color.orange;
                else if (lattice[i][j]<=0.998) c= new Color(255, 69, 0); // redOrange
 				else if (lattice[i][j]<=1.0) c=Color.red;
                else System.out.println ("Error here");
               /* if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                img.setRGB(i,j,c.getRGB());   
			}
		repaint();
		counter++;
	}
    

	
	public  BufferedImage getBM ()
	{
		BufferedImage result = new BufferedImage (ca.size,ca.size,BufferedImage.TYPE_INT_RGB);
		float[][] lattice = ca.getBM();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.1) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.3) c=Color.cyan;
                else if (lattice[i][j]<=0.4) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.7) c=Color.yellow;
				else if (lattice[i][j]<=0.8) c=Color.orange;  
				else if (lattice[i][j]<=0.9) c= new Color(255, 69, 0); // redOrange
				else if (lattice[i][j]<=1.0) c=Color.red;  
               /* if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                result.setRGB(i,j,c.getRGB());   
			}
		return result;

	}
	public void paintBM()
	{
		float[][] lattice = ca.getBM();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
                Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.1) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.3) c=Color.cyan;
                else if (lattice[i][j]<=0.4) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.7) c=Color.yellow;
				else if (lattice[i][j]<=0.8) c=Color.orange;  
				else if (lattice[i][j]<=0.9) c= new Color(255, 69, 0); // redOrange
				else if (lattice[i][j]<=1.0) c=Color.red;  
              /*  if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                img.setRGB(i,j,c.getRGB());   
			}
		repaint();
		counter++;
	}

	public BufferedImage getMDE ()
	{
		BufferedImage result = new BufferedImage (ca.size,ca.size,BufferedImage.TYPE_INT_RGB);
		float[][] lattice = ca.getMDE();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.15) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.25) c=Color.cyan;
                else if (lattice[i][j]<=0.3) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.8) c=Color.yellow;
				else if (lattice[i][j]<=0.995) c=Color.orange;
                else if (lattice[i][j]<=0.998) c= new Color(255, 69, 0); // redOrange
 				else if (lattice[i][j]<=1.0) c=Color.red;
               /* if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                result.setRGB(i,j,c.getRGB());   
			}
		return result;
	}
    
    public void paintMDE ()
	{
		float[][] lattice = ca.getMDE();
		for (int i=0;i<ca.size;i++)
			for (int j=0;j<ca.size;j++) {
				Color c = new Color (Color.blue.getRGB());
				if (lattice[i][j]<=0.01f) c=Color.white;
                else if (lattice[i][j]<=0.15) c= new Color(0, 0, 128); // darkBlue
                else if (lattice[i][j]<=0.2) c=Color.blue;
				else if (lattice[i][j]<=0.25) c=Color.cyan;
                else if (lattice[i][j]<=0.3) c= new Color(0, 128, 128); // darkCyan 
				else if (lattice[i][j]<=0.5) c=Color.green;
                else if (lattice[i][j]<=0.6) c= new Color(170, 255, 0); // limeGreen
				else if (lattice[i][j]<=0.8) c=Color.yellow;
				else if (lattice[i][j]<=0.995) c=Color.orange;
                else if (lattice[i][j]<=0.998) c= new Color(255, 69, 0); // redOrange
 				else if (lattice[i][j]<=1.0) c=Color.red;
                else System.out.println ("Error here");
              /*  if ((i==999) || (j==999) || (i==0) || (j==0) || (i==499) || (j==499) || (i==500) || (j==500) || (i==501) || (j==501) ) c=Color.white; // border  */
                img.setRGB(i,j,c.getRGB());   
			}
		repaint();
		counter++;
	}

	
	public void actionPerformed (ActionEvent e) {
		String res = e.getActionCommand();
		if ("play".equals(e.getActionCommand())) {
			if (play) {
				play=false;
				button.setText ("Play");
			}
			else {
				play=true;
				button.setText ("Pause");
			}
		}
		else if (res.compareTo("Movie")==0) {
			if (movie==false) movie=true;
			else movie=false;
		}
		else if (res.compareTo("openButton")==0) {
			fc.setFileFilter(new javax.swing.filechooser.FileFilter() {
					public boolean accept (File f) {
					return f.getName().toLowerCase().endsWith(".txt") || f.isDirectory();
					}
					public String getDescription() {
						return "prosCA config txt files";
						}

					});
            		int returnVal = fc.showOpenDialog(this);
            		if (returnVal == JFileChooser.APPROVE_OPTION) {
	                	File file = fc.getSelectedFile();
				try{
					label.setText ("Loading "+file.getCanonicalPath());
					//ca.readConfig(file.getCanonicalPath());
					ca.readConfig(file);
					ca.reset();
					label.setText (file.getCanonicalPath()+" loaded");
					rootDir=file.getParent();
				} catch (Exception ex) {};
				//try{System.out.println (file.getCanonicalPath());} catch (Exception ex) {};
			}
        	        //This is where a real application would open the file.
        	} else { 
			defaultView =views.getSelectedIndex(); 
		}
		repaint();
	}
	public void mouseExited (MouseEvent e) {}	
	public void mouseEntered (MouseEvent e) {}
	public void mouseReleased (MouseEvent e) {}
	public void mousePressed (MouseEvent e) {}


	public static void main(String args[]) {
		int mag=1;
		int maxTS=-1;
		boolean movie;
	
		System.out.println ("# Vis version:"+Vis3.version);
		System.out.println ("# CA version:"+prosCA3.version);
		
		// Let's deal now with the main window
		JFrame f = new JFrame("Prostate CA Visualisation "+prosCA3.version);
	        f.addWindowListener(new WindowAdapter() {
	        	public void windowClosing(WindowEvent e) {System.exit(0);}
	        });
		JLabel label = new JLabel ("Timesteps");

		// Now the views list
		String[] theViews = {"Cells","TGF-B","MDE","BM"};
		JComboBox views = new JComboBox(theViews);
		views.setSelectedIndex(0);
		
		
		JButton openButton = new JButton("Set a parameter file...");
		JButton play = new JButton ("Play");
		JCheckBox cb = new JCheckBox ("Movie",false);
		Vis3 m = new Vis3 (views,play,label,cb);
		views.addActionListener(m);
	        openButton.addActionListener(m);
		play.setActionCommand("play");
		openButton.setActionCommand ("openButton");
		play.addActionListener (m);
		cb.addActionListener(m);
		JPanel buttonPanel = new JPanel(); //use FlowLayout
		buttonPanel.add(cb);
		buttonPanel.add(openButton);
		buttonPanel.add(play);
		buttonPanel.add(views);
	
        f.getContentPane().add("Center", m);
		f.getContentPane().add("North",label);
		//f.getContentPane().add("South",views);
		f.getContentPane().add("South",buttonPanel);
		
		
	        f.pack();
	        f.setSize(new Dimension(2*m.ca.mag*m.size,2*m.ca.mag*m.size));
	        f.show();
		
		boolean finished=false;
		int ts=0;
		while (ts< m.ca.timesteps) {
            System.out.print ("");
			if (m.play) {
				m.nextTimeStep();
				label.setText ("Timestep: "+ts);
				ts++;
			}
		} 
	    }

}
