/*
 * Algorithm developed by: Ziv Frankenstein, Ph.D. 
 */


import java.io.*;
import java.util.*;

 
public class model {
	public static String version="30 april 2012";
	MersenneTwisterFast random;
	public int mag=1;
	public int timesteps=-1;
	public float [][] TGFB=null;
	public float [][] MDE=null; 
	public float [][] BM=null; 
	public float [][] TGFBProduction=null;
	public float [][] MDEProduction=null;
	int [][] Cells=null; 
    int [][] outBasLum=null; 
	arrayFF cellList=null; 
	int [][] EpithelialCells=null; 
	float [][] Age=null; 
	public int numCells=0;
	boolean finishedRun=false;
	public int totalTumourCells=0; 
	
	int size = 1000; /
	int timestep=0; 
	float[] migrationProbabilityCoefficient = {0.0f, 0.0f, 0.0f, 0.0001f, 0.0f, 0.1f }; 
	float[] proliferatingProbabilityCoefficient = { 0.0f, 0.0f, 0.0f, 0.5f, 0.0f, 0.0f }; 
	float[] thresholdTGFB = { 
	 0.0f,	 /* no cells */ 
	 0.5f, /* basal cells */
	 0.5f, /* normal luminal */
	 0.0f, /* normal stromal  */
	 0.5f  
	};
	float[] thresholdMDE = { 
	 0.0f, /* no cells */
	 0.5f, /* basal cells */
	 0.5f, /* normal luminal */
	 0.5f, /* normal stromal */ 
	 0.5f  
	};
	float tgfbDecayRate = 0.001f;
	float probabilityOfMonocyteAppearing = 0.0f; 
	float probabilityOfApoptosisBasal = 0.01f;
	float probabilityOfNecrosis = 0.1f;
	float abnormalLuminalCellMDEProductionRateConstant = 0.1f; 
	public float abnormalLuminalCellTGFBProductionRateConstant = 1.0f;
	float basalCellBMProductionRateConstant = 0.1f;
	float kDe = 0.03f;
 	float mdeDecayRate = 0.001f;
	float bmDecayRate = 0.1f;
	float basalCellTGFBproductionRate = 0.01f;
	float luminalCellTGFBconsumptionRateConstant = 0.5f;
	float kDt= 0.02f;
	float h=1.0f;
	public float wiggleRoom=0.1f;
    public float avgTGFBProduction;
    public float avgMDEProduction;
    public float OutDuctTGF=-1; 
    public float OutDuctMDE=-1; 
    public int OutDuctCoorX=-1; 
    public int OutDuctCoorY=-1; 
    public int gapX=100000; 
    public int gapY=100000; 

	int xCenter = size / 2 - 1;
	 int yCenter = size / 2 - 1;


	// Model Parameters
	float initTGFB=0.0f;
	float initMDE=0.00f; 
	float initBM=1.0f; 
	int mutTS=10; 
	float stromalCellTGFBproductionRate = 0.02f;
	float stromalMembraneProductionRate = 0.1f;	
    float stromalTGFBConsumptionRate = 0.1f;
	float densityMonocytes=0.01f;
	float densityStroma=0.1f;
	float densityECM=0.1f;
	float stromalTGFBswitch = 0.01f;
	float ECMTGFBconsumption = 0.0001f;
	float abnormalLuminalCellApoptosisThreshold=0.1f;
	float probBasalCellMitosis=0.1f;
    float mutationRate = 0.0f;  
	public String OutName=""; 

        
	public model ()
	{
		int time = (int) System.currentTimeMillis();
		random = new MersenneTwisterFast (time);
		readConfig(new File("parameters.txt"));
        DataPrint();  
		reset();
	}
	
	public void reset ()
	{
		resetTGFB();
		resetMDE();
		resetCells("ic1000.txt");
    }
    
    void readConfig (File file)
    {
            try {
                String filename = file.getCanonicalPath();
                System.out.println ("Reading: "+filename);
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String str;
		        while ((str = in.readLine()) != null) {

                    if ((str.length()>0) && (str.charAt(0)!='#')){
                        StringTokenizer st = new StringTokenizer (str);
                        String name = st.nextToken();
                        st.nextToken();
                        String value = st.nextToken();
                        if (name.compareToIgnoreCase("size")==0) size=Integer.parseInt(value);
                            else if (name.compareToIgnoreCase("initTGFB")==0) initTGFB= Float.parseFloat (value);
                                else if (name.compareToIgnoreCase("initMDE")==0) initMDE= Float.parseFloat (value);
                                    else if (name.compareToIgnoreCase("initBM")==0) initBM=Float.parseFloat (value);
                                        else if (name.compareToIgnoreCase("stromalTGFBproduction")==0) stromalCellTGFBproductionRate =Float.parseFloat (value);
                                            else if (name.compareToIgnoreCase("tgfbDecayRate")==0) tgfbDecayRate =Float.parseFloat (value);
                                                else if (name.compareToIgnoreCase("mdeDecayRate")==0) mdeDecayRate=Float.parseFloat (value); 
                                                    else if (name.compareToIgnoreCase("bmDecayRate")==0) bmDecayRate=Float.parseFloat (value);
                                                        else if (name.compareToIgnoreCase("kDe")==0) kDe=Float.parseFloat (value);
                                                            else if (name.compareToIgnoreCase("kDt")==0) kDt=Float.parseFloat (value);
                                                                else if (name.compareToIgnoreCase("probabilityOfApoptosisBasal")==0) probabilityOfApoptosisBasal=Float.parseFloat (value);
                                                                    else if (name.compareToIgnoreCase("probabilityOfNecrosis")==0) probabilityOfNecrosis=Float.parseFloat (value);
                                                                        else if (name.compareToIgnoreCase("abnormalLuminalCellMDEProductionRateConstant")==0) abnormalLuminalCellMDEProductionRateConstant=Float.parseFloat (value);
                                                                            else if (name.compareToIgnoreCase("basalCellBMProductionRateConstant")==0) basalCellBMProductionRateConstant=Float.parseFloat (value);
                                                                                else if (name.compareToIgnoreCase("basalCellTGFBproductionRate")==0) basalCellTGFBproductionRate=Float.parseFloat (value);
                                                                                    else if (name.compareToIgnoreCase("luminalCellTGFBconsumptionRateConstant")==0) luminalCellTGFBconsumptionRateConstant =Float.parseFloat (value);
                                                                                        else if (name.compareToIgnoreCase("abnormalLuminalCellTGFBProductionRateConstant")==0) abnormalLuminalCellTGFBProductionRateConstant =Float.parseFloat (value);
                                                                                            else if (name.compareToIgnoreCase("stromalMembraneProductionRate")==0) stromalMembraneProductionRate = Float.parseFloat (value);
                                                                                                else if (name.compareToIgnoreCase("stromalTGFBConsumptionRate")==0) stromalTGFBConsumptionRate = Float.parseFloat (value);
                                                                                                    else if (name.compareToIgnoreCase("densityMonocytes")==0) densityMonocytes = Float.parseFloat (value);
                                                                                                        else if (name.compareToIgnoreCase("densityStroma")==0) densityStroma = Float.parseFloat (value);
                                                                                                            else if (name.compareToIgnoreCase("densityECM")==0) densityECM = Float.parseFloat (value);
                                                                                                                else if (name.compareToIgnoreCase("stromalTGFBswitch")==0) stromalTGFBswitch = Float.parseFloat (value);
                                                                                                                    else if (name.compareToIgnoreCase("magnification")==0) mag = Integer.parseInt (value);
                                                                                                                        else if (name.compareToIgnoreCase("ECMTGFBconsumption")==0) ECMTGFBconsumption = Float.parseFloat (value);
                                                                                                                            else if (name.compareToIgnoreCase("abnormalLuminalCellApoptosisThreshold")==0) abnormalLuminalCellApoptosisThreshold = Float.parseFloat (value);
                                                                                                                                else if (name.compareToIgnoreCase("probBasalCellMitosis")==0) probBasalCellMitosis = Float.parseFloat (value);
                                                                                                                                    else if (name.compareToIgnoreCase("timesteps") ==0) timesteps = Integer.parseInt (value);
                                                                                                                                        else if (name.compareToIgnoreCase("wiggleRoom") ==0) wiggleRoom = Float.parseFloat (value);
                                                                                                                                           else if (name.compareToIgnoreCase("OutName") ==0) OutName = value;  
                                                                                                                                            else if (name.compareToIgnoreCase("mutationRate") ==0) mutationRate = Float.parseFloat (value);  
                                                                                                                                              else System.err.println (name+": Unrecognised field");
                                                                                                                                                }
        		}
        		in.close();
            } catch (Exception e) {
                System.err.println ("Configuration file not found");
                e.printStackTrace();
            }
        }

        public void DataPrint()                    

        {          
            try{
                Writer outputData = null;
                
                String strDirectoy = OutName;   
                boolean success = (new File(strDirectoy)).mkdir(); 

                File file = new File(OutName+"/"+OutName+"Data.out");
                outputData = new BufferedWriter(new FileWriter(file));

                FileInputStream fstream = new FileInputStream("parameters.txt");
                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                while ((strLine = br.readLine()) != null)   {
                    outputData.write(strLine);
                    outputData.write("\n");
                }


                in.close();  
                outputData.close();              

            }catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }
        }
                
    
	final int distance (int x, int y, int i, int j)
	{
		double dis=Math.sqrt((x-i)*(x-i)+(y-j)*(y-j));

		return (int)Math.round(dis);
	}

	final int[] convertCoordinates(int x, int y)
       {

               if (x < 0) x = size - 1;
               else if (x > size - 1) x = 0;
               if (y < 0) y = size - 1;
               else if (y > size - 1) y = 0;
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }
           
        void cancerInitiation()
                   {
            boolean done=false;
            while (done==false) {
                int x = 496;
                int y = 494;

                if ( (Cells[x][y]==2) ) { 
                    done=true;
                    Cells[x][y]=3;
                    TGFBProduction[x][y]=abnormalLuminalCellTGFBProductionRateConstant;
                    MDEProduction[x][y]=abnormalLuminalCellMDEProductionRateConstant;
                    maxTGFBProd=TGFBProduction[x][y];
                    maxMDEProd=MDEProduction[x][y];
                }
            }
        }
    
    
	public void nextTimeStep ()
    
	{
		totalBirths=0;
        totalDeaths=0;
                
		if (timestep==mutTS) cancerInitiation();
		for (int j=0;j<5;j++) {
			iterateTGFB();
			iterateMDE();
			iterateBM();
		}
		iterateCells();
        
		timestep++;
        
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Age[i][j]+=0.1f;
                totalMembrane+=BM[i][j];
                
				if ((Cells[i][j]==1) || (Cells[i][j]==2)) totalEpithelialCells++;
				if (Cells[i][j]==3) {
					totalTumourCells++;
					totalTGFBProduction+=TGFBProduction[i][j];
                    totalMDEProduction+=MDEProduction[i][j];
				}

                if (  (Cells[i][j]==3) && (  TGFBProduction[i][j] > 0.01 )) {
                    SitetotalTumourCellsTGF++;
					siteTGFBProduction+=TGFBProduction[i][j];
                    
                }
                if (  (Cells[i][j]==3) && (  MDEProduction[i][j] > 0.01 )) {
                    SitetotalTumourCellsMDE++;
					siteMDEProduction+=MDEProduction[i][j];
                    
                }
                if ((OutDuctTGF==-1) && (outBasLum[i][j] !=-333) && (  TGFB[i][j] > 0.01 )) {
                    outBasLum[i][j] = 111;
                }
                if ( (outBasLum[i][j] != 111) && (  TGFB[i][j] > 0.01 )){
                    SitetotalTgfCellsGen++;
					SiteTGFBGen+=TGFB[i][j];
                    
                }
                if (  (  TGFB[i][j] > 0.01 )){
                    SitetotalTgfCellsGenAll++;
					SiteTGFBGenAll+=TGFB[i][j];
                }

                if (  (  MDE[i][j] > 0.01 )){
                    SitetotalMdeCellsGen++;
					SiteMDEGen+=MDE[i][j];
                }
                    
                if ( (Cells[i][j]==6) ) {
                    activeStromaWhereTumour++;
                }
                if ( (Cells[i][j]==5) ) {
                    monocyteWhereTumour++;
                }
                if ( (Cells[i][j]==4) ) {
                    restStromaWhereTumour++;
                }
                
			}

        
        SiteavgTGFBProduction=siteTGFBProduction/(SitetotalTumourCellsTGF);
		SiteavgMDEProduction=siteMDEProduction/SitetotalTumourCellsMDE;
        SiteavgTGFBGen=SiteTGFBGen/(SitetotalTgfCellsGen);
        SiteavgTGFBGenAll=SiteTGFBGenAll/(SitetotalTgfCellsGenAll);
        SiteavgMDEGen=SiteMDEGen/(SitetotalMdeCellsGen);
        
		
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
                
                if (  (Cells[i][j]==3) && (  TGFBProduction[i][j] > 0.01 )) {
                    SitestdDevTGFBprod+=(SiteavgTGFBProduction-TGFBProduction[i][j])*(SiteavgTGFBProduction-TGFBProduction[i][j]);
                }    
                if (  (Cells[i][j]==3) && (  MDEProduction[i][j] > 0.01) ) {
					SitestdDevMDEprod+=(SiteavgMDEProduction-MDEProduction[i][j])*(SiteavgMDEProduction-MDEProduction[i][j]);
                }
                if ((  TGFB[i][j] > 0.01 ) && (outBasLum[i][j] != 111) ) {
                    SitestdDevTGFBGen+=(SiteavgTGFBGen-TGFB[i][j])*(SiteavgTGFBGen-TGFB[i][j]);
                } 
                if (  TGFB[i][j] > 0.01 ) {
                    SitestdDevTGFBGenAll+=(SiteavgTGFBGenAll-TGFB[i][j])*(SiteavgTGFBGenAll-TGFB[i][j]);
                }    
                if (  MDE[i][j] > 0.01 ) {
                    SitestdDevMDEGen+=(SiteavgMDEGen-MDE[i][j])*(SiteavgMDEGen-MDE[i][j]);                
                }
            }
        
        SitestdDevTGFBprod=(float)Math.sqrt(SitestdDevTGFBprod/SitetotalTumourCellsTGF);
        SitestdDevMDEprod=(float)Math.sqrt((float)(SitestdDevMDEprod/SitetotalTumourCellsMDE));
        SitestdDevTGFBGen=(float)Math.sqrt(SitestdDevTGFBGen/SitetotalTgfCellsGen);
        SitestdDevTGFBGenAll=(float)Math.sqrt(SitestdDevTGFBGenAll/SitetotalTgfCellsGenAll);
        SitestdDevMDEGen=(float)Math.sqrt((float)(SitestdDevMDEGen/SitetotalMdeCellsGen));
        
        if (totalTumourCells > 0) {    
            
            for (int i=0;i<size;i++)
                for (int j=0;j<size;j++) {
                    if ( ((size-i) <= (i-0)) && ((size-i) < gapX) && (Cells[i][j]==3) ) {
                        gapX = size-i;
                        OutDuctCoorX = i;
                    }
                    if ( ((size-i) >= (i-0)) && ((i-0) < gapX) && (Cells[i][j]==3) ) {
                        gapX = i-0;
                        OutDuctCoorX = i;
                    }
                    if ( ((size-j) <= (j-0)) && ((size-j) < gapY) && (Cells[i][j]==3) ) {
                        gapY = size-j;
                        OutDuctCoorY = j;
                    }
                    if ( ((size-j) >= (j-0)) && ((j-0) < gapY) && (Cells[i][j]==3) ) {
                        gapY = j-0;
                        OutDuctCoorY = j;
                    }
                }
            
            if (gapX >= gapY){
                for (int i=0;i<size;i++)
                    if ( (Cells[i][OutDuctCoorY]==3) ) OutDuctCoorX = i;
            }
            if (gapX <= gapY){
                for (int j=0;j<size;j++)
                    if ( (Cells[OutDuctCoorX][j]==3) ) OutDuctCoorY = j;
            }
            if ( (outBasLum[OutDuctCoorX][OutDuctCoorY]==-333) && (Cells[OutDuctCoorX][OutDuctCoorY]==3) ) {
                OutDuctTGF=TGFBProduction[OutDuctCoorX][OutDuctCoorY]; 
                OutDuctMDE=MDEProduction[OutDuctCoorX][OutDuctCoorY]; 
            }
         }   
        
        
        BufferedReader reader = null;          
        BufferedWriter writer = null;
        ArrayList list = new ArrayList();
        
        try {
            
            reader = new BufferedReader(new FileReader(OutName+"/"+OutName+"Data.out"));
            String tmp;
            while ((tmp = reader.readLine()) != null)
                list.add(tmp);
            reader.close();
            
            if (timestep==1) list.add ("timestep\ttotalEpithelial\ttotalTumourCells\tSitetotalTumourCells\tSitetotalTumourCellsDucts\tSiteavgTGFBProduction\tSiteavgMDEProduction\tSiteavgTGFB\tSiteavgMDE\tSiteavgTGFBProductionDucts\tSiteavgMDEProductionDucts\tSiteavgTGFBDucts\tSiteavgMDEDucts\tSitetotalTumourCellsGen\tSiteavgTGFBProductionGen\tSiteavgMDEProductionGen\tSiteavgTGFBGen\tSiteavgMDEGen\tSitestdDevTGFBprod\tSitestdDevMDEprod\tSitestdDevTGFBGen\tSitestdDevMDEGen\tOutDuctTGF\tOutDuctMDE\tOutDuctCoorX\tOutDuctCoorY\tactiveStromaWhereTumour\tmonocyteWhereTumour\trestStromaWhereTumour\tMDE\tMembrane\tTGFBoutDuct\tTGFB\tTGFBProduction\tstdDe\tMDEProduction\tBirths\tDeaths\tmaxTGFBProd\tminTGFBProd\tSiteavgTGFBGenAll\tSitestdDevTGFBGenAll\n");
            
            list.add (timestep+"\t"+totalEpithelialCells+"\t"+totalTumourCells+ "\t"+SitetotalTumourCells+ "\t"+SitetotalTumourCellsDucts+ "\t"+ SiteavgTGFBProduction+ "\t"+SiteavgMDEProduction+ "\t"+SiteavgTGFB+ "\t"+SiteavgMDE+ "\t"+SiteavgTGFBProductionDucts+ "\t"+SiteavgMDEProductionDucts+ "\t"+SiteavgTGFBDucts+"\t"+SiteavgMDEDucts+ "\t"+SitetotalTumourCellsGen+ "\t"+SiteavgTGFBProductionGen+ "\t"+SiteavgMDEProductionGen+ "\t"+SiteavgTGFBGen+ "\t"+SiteavgMDEGen+ "\t"+SitestdDevTGFBprod+ "\t"+SitestdDevMDEprod+ "\t"+SitestdDevTGFBGen+ "\t"+SitestdDevMDEGen +"\t"+OutDuctTGF+ "\t"+OutDuctMDE+"\t"+OutDuctCoorX+"\t"+OutDuctCoorY +"\t"+activeStromaWhereTumour+"\t"+monocyteWhereTumour+"\t"+restStromaWhereTumour+"\t"+SiteMDEGen+"\t"+totalMembrane+"\t"+SiteTGFBGen+"\t"+SiteTGFBGenAll+"\t"+siteTGFBProduction+"\t"+stdDevTGFB+"\t"+siteMDEProduction+"\t"+totalBirths+"\t"+totalDeaths+"\t"+maxTGFBProd+"\t"+minTGFBProd+"\t"+SiteavgTGFBGenAll+"\t"+SitestdDevTGFBGenAll);
            
            writer = new BufferedWriter(new FileWriter(OutName+"/"+OutName+"Data.out"));
            for (int i = 0; i < list.size(); i++)
                writer.write(list.get(i) + "\r\n");
            writer.close();
            if ((OutDuctCoorX==0) || (OutDuctCoorX==1000) || (OutDuctCoorY==0) || (OutDuctCoorY==1000)) System.exit(0); // kill when tumor get into grid border
            
		} catch (IOException e) {
			e.printStackTrace();
        } 
    }

	void resetTGFB ()
	{
		TGFB= new float [size][size];
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				TGFB[i][j]=initTGFB;
			}

	}
	
	void resetMDE()
	{
		MDE= new float [size][size];
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				MDE[i][j]=initMDE;
			}

	}
	
                   
    void fillBM (int x, int y)
    {
        int[] p = new int [2];
        
        p=convertCoordinates (x+1,y-1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x+1,y);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x+1,y+1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x,y-1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x,y+1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x-1,y-1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x-1,y);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
        p=convertCoordinates (x-1,y+1);
        if ( (Cells[p[0]][p[1]]<3) && (p[0]>0) && (p[0]<999) && (p[1]>0) && (p[1]<999) ) BM[p[0]][p[1]]=initBM;
    }
                   
    void fillLuminal (int x, int y)
    {
            int[] p = new int [2];
            
            p=convertCoordinates (x+1,y-1);
            if (Cells[p[0]][p[1]]==0) Cells[p[0]][p[1]]=2;
                p=convertCoordinates (x+1,y);
                if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                    p=convertCoordinates (x+1,y+1);
                    if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                        p=convertCoordinates (x,y-1);
                        if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                            p=convertCoordinates (x,y+1);
                            if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                                p=convertCoordinates (x-1,y-1);
                                if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                                    p=convertCoordinates (x-1,y);
                                    if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
                                        p=convertCoordinates (x-1,y+1);
                                        if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=2;
    }
                   
    void fillBasal (int x, int y)
    {
            int[] p = new int [2];
            
            p=convertCoordinates (x+1,y-1);
            if (Cells[p[0]][p[1]]==0) Cells[p[0]][p[1]]=1;
                p=convertCoordinates (x+1,y);
                if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                    p=convertCoordinates (x+1,y+1);
                    if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                        p=convertCoordinates (x,y-1);
                        if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                            p=convertCoordinates (x,y+1);
                            if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                                p=convertCoordinates (x-1,y-1);
                                if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                                    p=convertCoordinates (x-1,y);
                                    if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
                                        p=convertCoordinates (x-1,y+1);
                                        if (Cells[p[0]][p[1]]==0)  Cells[p[0]][p[1]]=1;
    }
                   
    void resetCells (String filename)
                          {
                          try {
                          BufferedReader in = new BufferedReader(new FileReader(filename));
                          String str;
                          str=in.readLine();
                          StringTokenizer stSize = new StringTokenizer (str);
                          int size = Integer.parseInt (stSize.nextToken());
                          int centre = size/2;
                          Cells = new int [size][size];
                          EpithelialCells = new int [size][size];
                          outBasLum = new int [size][size];
                          Age = new float [size][size];
                          BM=new float [size][size];
                          MDEProduction = new float [size][size];
                          TGFBProduction = new float[size][size];
                          for (int i=0;i<size;i++)
                          for (int j=0;j<size;j++) {
                          Cells[i][j]=0;
                          BM[i][j]=densityECM;
                          MDEProduction[i][j]=0;
                          TGFBProduction[i][j]=0;
                          
                          }
                          
                          for (int j=0;j<size;j++) {
                          str=in.readLine();
                          StringTokenizer st = new StringTokenizer (str);
                          for (int i=0;i<size;i++) Cells[i][j] = Integer.parseInt (st.nextToken());
                          }
                          
                          for (int i=0;i<size;i++) 
                          for (int j=0;j<size;j++) 
                          if (Cells[i][j]==7) fillLuminal (i,j);
                          
                          for (int i=0;i<size;i++)
                          for (int j=0;j<size;j++) 
                          if (Cells[i][j]==2) fillBasal (i,j);
                          
                          for (int i=0;i<size;i++)
                          for (int j=0;j<size;j++) {
                          if ((Cells[i][j]==1) || (Cells[i][j]==2)) EpithelialCells[i][j]=Cells[i][j];
                          Age[i][j] = random.nextFloat();
                          }
                          
                          for (int i=0;i<size;i++)
                          for (int j=0;j<size;j++) {
                          if ((Cells[i][j]==1)|| (Cells[i][j]==2)){
                          fillBM (i,j);
                          if (Cells[i-1][j]<3) fillBM (i-1,j);
                          if (Cells[i+1][j]<3) fillBM(i+1,j);
                          if (Cells[i][j-1]<3) fillBM(i,j-1);
                          if (Cells[i][j+1]<3) fillBM (i,j+1);
                          if (Cells[i-1][j-1]<3) fillBM (i-1,j-1);
                          if (Cells[i-1][j+1]<3) fillBM (i-1,j+1);
                          if (Cells[i+1][j-1]<3) fillBM (i+1,j-1);
                          if (Cells[i+1][j+1]<3) fillBM (i+1,j+1);
                          }
                          }
                          
                          int numMonocytes=(int) (size*size*densityMonocytes); 
                          while (numMonocytes!=0) {
                          int posX=random.nextInt(size);
                          int posY=random.nextInt(size);
                          if ( (Cells[posX][posY]==0) && (BM[posX][posY]!=initBM) ) {
                          Cells[posX][posY]=5;
                          MDEProduction[posX][posY]=0; 
                          TGFBProduction[posX][posY]=0; 
                          Age[posX][posY]=0;
                          numMonocytes--;
                          outBasLum[posX][posY]=-333; 
                          }	
                          }
                          int numStroma=(int) (size*size*densityStroma);  
                          while (numStroma!=0) {
                          int posX=random.nextInt(size);
                          int posY=random.nextInt(size);
                          if ( (Cells[posX][posY]==0) && (BM[posX][posY]!=initBM) ) {
                          Cells[posX][posY]=4;
                          MDEProduction[posX][posY]=0; 
                          TGFBProduction[posX][posY]=0; 
                          Age[posX][posY]=0;
                          numStroma--;
                          outBasLum[posX][posY]=-333; 
                          }	
                          }
                          
                          for (int i=0;i<size;i++)
                          for (int j=0;j<size;j++) {
                          if (Cells[i][j]==0) outBasLum[i][j]=-333;  
                          if (Cells[i][j]==7) BM[i][j]=0;
                          if (Cells[i][j]==7) Cells[i][j]=0;
                          }
                          
                          } catch (IOException e) {
                          e.printStackTrace();
                          System.exit (-1);
                          }
                          }
                   
   

	public boolean iterateCells()
 	{
		if (cellList==null) cellList = new arrayFF (size*size);
        int inTumorCells = 0;
        
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++) {
				if (Cells[i][j] > 0)  { 
					int[] p = new int[2];
					p[0] = i; p[1] = j;
					cellList.add(p);
				}
            }
                
		while (cellList.size() != 0) {
			int randomElemIndex=0;
			if (cellList.size()>1) randomElemIndex = random.nextInt(cellList.size()-1);
	   		int[] point = (int[])cellList.get(randomElemIndex); 
		   	int rI = point[0];
		   	int rJ = point[1];
            
			cellList.remove(randomElemIndex); 
		   	int cell = Cells[rI][rJ];
            
			int [] newP = getNewCellLoc(cell,rI,rJ);
            cell = Cells[rI][rJ];
            
			if (newP[0] == rI && newP[1] == rJ) {
				boolean didProliferate = isProliferating(cell,rI,rJ); 
				if (didProliferate==false) {
					switch (Cells[rI][rJ]) {
						case 1:
                            if (random.nextFloat()<=probabilityOfApoptosisBasal*MDE[rI][rJ]) {
                                Cells[rI][rJ]=0;
                                MDEProduction[rI][rJ]=0; 
                                TGFBProduction[rI][rJ]=0; 
                            }
							break;
						case 2:
                            if (random.nextFloat()<=getCancerNeighbours(rI,rJ)*probabilityOfNecrosis) {
                                Cells[rI][rJ]=0;
                                MDEProduction[rI][rJ]=0; 
                                TGFBProduction[rI][rJ]=0; 
                            }	
							break;
						case 3: 
							Age[rI][rJ] += getStromaMonocyteNeighbours(rI, rJ) * TGFB[rI][rJ]; 
							if (TGFB[rI][rJ]<= abnormalLuminalCellApoptosisThreshold)
								if (random.nextFloat()<=(abnormalLuminalCellApoptosisThreshold-TGFB[rI][rJ])/abnormalLuminalCellApoptosisThreshold) {
									Cells[rI][rJ]=0;
                                    MDEProduction[rI][rJ]=0; 
                                    TGFBProduction[rI][rJ]=0; 
									totalDeaths++;
									System.err.println ("Tumour cell dies due to insufficient TGFB: "+TGFB[rI][rJ]);
								}
							break;
                    }
				}
		   	} else {
				int x, y;
				if (newP[0]==size) x = 0;
                else if (newP[0]==-1) x=size-1;
                else  x=newP[0];
				if (newP[1]==size) y = 0;
                else if (newP[1]==-1) y=size-1;
                else  y=newP[1];
                if ( ((x==0) || (x==1000) || (y==0) || (y==1000)) && (cell==3) ) {
                    OutDuctTGF=TGFBProduction[rI][rJ]; 
                    OutDuctMDE=MDEProduction[rI][rJ]; 
                    OutDuctCoorX=x;
                    OutDuctCoorY=y;
                }
                if ( (outBasLum[x][y]==-333) && (cell==3) && (OutDuctTGF==-1)) {
                    OutDuctTGF=TGFBProduction[rI][rJ]; 
                    OutDuctMDE=MDEProduction[rI][rJ]; 
                    OutDuctCoorX=x;
                    OutDuctCoorY=y;
                }
                Cells[x][y] = cell;
				Age[x][y]=Age[rI][rJ];
                MDEProduction[x][y]=0;
                TGFBProduction[x][y]=0;
                if (cell==3) {  
                    MDEProduction[x][y]=MDEProduction[rI][rJ]; 
                if ((cell ==1) || (cell ==3) || (cell ==4)) {   
                    TGFBProduction[x][y]=TGFBProduction[rI][rJ]; 
                }   
				Cells[rI][rJ] = 0;
				Age[rI][rJ]=0;
                MDEProduction[rI][rJ]=0; 
                TGFBProduction[rI][rJ]=0; 
	   		}
            
	  	}
        return true;
    }

public int getStromaMonocyteNeighbours(int  rI, int rJ)
{
	int stromaSites = 0;
	for (int dI = -1; dI <= 1; dI++)
		if (((rI + dI) >= 0) && ((rI + dI) < size))
			for (int dJ = -1; dJ <= 1; dJ++) 
				if ((rJ + dJ >= 0) && ((rJ + dJ) < size)) 
					if ((Cells[rI + dI][rJ + dJ] == 4) || (Cells[rI + dI][rJ + dJ] == 5)) stromaSites++; 
	return stromaSites;

}


public int[] getNewCellLoc(int cell, int rI, int rJ)
 {

	int[] newLoc = { rI, rJ }; 

	if ((cell == 3 && BM[rI][rJ] < 0.01) || (cell == 5))
	{
		
		if (random.nextFloat() <= migrationProbabilityCoefficient[cell]) {
			switch (cell) {
				case 3:  
					return determineLoc(cell, rI, rJ,false);
				case 5: 
					return determineLoc(cell, rI, rJ,true);
			}
		} 
	} 
     
     return newLoc; 
}

public int getCancerNeighbours (int rI, int rJ) {
	int cancerSites = 0;
  	for (int dI = -1; dI <= 1; dI++)
   		if (((rI + dI) >= 0) || ((rI + dI) <= size - 1)) 
    			for (int dJ = -1; dJ <= 1; dJ++)
     				if ((rJ + dJ >= 0) || ((rJ + dJ) <= size - 1)) 
						if (Cells[rI + dI][rJ + dJ] == 3) cancerSites++; 
	return cancerSites;
}

public int[] determineLoc(int cell, int rI, int rJ, boolean chemotaxisTGFB) {
	LinkedList vacantSites = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	if (rI+1==size) {tp1[0] = 0; tp1[1] = rJ;}
		else {tp1[0] = rI + 1; tp1[1] = rJ;}
	if (Cells[tp1[0]][tp1[1]]==0) vacantSites.add(tp1);
	if (rI-1==-1) {tp2[0] = size-1 ; tp2[1] = rJ;}
		else {tp2[0] = rI -1 ; tp2[1] = rJ;}
	if (Cells[tp2[0]][tp2[1]]==0) vacantSites.add(tp2);
	if (rJ+1==size) {tp3[0] = rI; tp3[1] = 0;}
		else {tp3[0] = rI; tp3[1] = rJ+1;}
	if (Cells[tp3[0]][tp3[1]]==0) vacantSites.add(tp3);
	if (rJ-1==-1) {tp4[0] = rI; tp4[1] = size-1;}
		else {tp4[0] = rI; tp4[1] = rJ-1;}
	if (Cells[tp4[0]][tp4[1]]==0) vacantSites.add(tp4);
	
	if ((vacantSites.size() > 0) && chemotaxisTGFB) {
		int[] p = (int [])vacantSites.get(0);
		float k=0.01f;
		float h=0.005f;
		float D=0.00054f;
		float ro=0.01f;
		float P0, P1, P2, P3, P4;
		P0=P1=P2=P3=P4=0;
		
		float tgfbA, tgfbB, tgfbC, tgfbD;
		if (rI+1>=size) tgfbA = TGFB[0][rJ];
			else tgfbA=TGFB[rI+1][rJ];
		if (rI-1<0) tgfbB = TGFB[size-1][rJ];
			else tgfbB  = TGFB[rI-1][rJ];
		if (rJ+1>=size) tgfbC = TGFB[rI][0];
			else tgfbC=TGFB[rI][rJ+1];
		if (rJ-1<0) tgfbD = TGFB[rJ][size-1];
			else tgfbD  = TGFB[rI][rJ-1];
		P1 = (k*D/(h*h))-((k*ro)/(4*h*h))*(tgfbA-tgfbB);
		if (P1<0) P1=0;
			else if ((rI-1==-1) && (Cells[size-1][rJ]!=0)) P1=0;
			else if ((rI-1>=0) && (Cells[rI-1][rJ]!=0)) P1=0;
		P2 = (k*D/(h*h))+((k*ro)/(4*h*h))*(tgfbA-tgfbB);
		if (P2<0) P2=0;
			else if ((rI+1>=size) && (Cells[0][rJ]!=0)) P2=0;
			else if ((rI+1<size) && (Cells[rI+1][rJ]!=0)) P2=0;
		P3 = (k*D/(h*h))-((k*ro)/(4*h*h))*(tgfbC-tgfbD);
		if (P3<0) P3=0;
			else if ((rJ-1==-1) && (Cells[rI][size-1]!=0)) P3 =0;
			else if ((rJ-1>=0) && (Cells[rI][rJ-1]!=0)) P3=0;
		P4 = (k*D/(h*h))+((k*ro)/(4*h*h))*(tgfbC-tgfbD);
		if (P4<0) P4=0;
			else if ((rJ+1==size) && (Cells[rI][0]!=0)) P4=0;
			else if ((rJ+1<size) && (Cells[rI][rJ+1]!=0)) P4=0;
		P0 = 1-(4*k*D/(h*h))-(k*ro/(h*h))*(tgfbA+tgfbB-4*TGFB[rI][rJ]+tgfbC+tgfbD);
		if (P0<0) P0=0;
		float total = P0+P1+P2+P3+P4;
		if (total!=0){
			P0=P0/total;
			P1=P1/total;
			P2=P2/total;
			P3=P3/total;
			P4=P4/total;
		}
		float pos = random.nextFloat();
		int[] result = new int [2];
		if (pos<=P1) {
			if (rI-1==-1) {p[0] = size-1; p[1] = rJ;}
			 else {p[0] = rI-1; p[1] = rJ;} 
			 return p;
		}
		else if (pos<=P1+P2)  {
			if (rI+1==size) {p[0] = 0; p[1] = rJ;}
			else {p[0] = rI+1; p[1] = rJ;} 
            return p;
		}
		else if (pos<=P1+P2+P3)  {
			if (rJ-1==-1) {p[0] = rI; p[1] = size-1;}
			else {p[0] = rI; p[1] = rJ-1;} 
            return p;
		}
		else if (pos<=P1+P2+P3+P4)  {
			if (rJ+1==size) {p[0] = rI; p[1] = 0;} 
			else {p[0] = rI; p[1] = rJ+1;}
            return p;
		}
		else {
			p[0] = rI; p[1] = rJ;
			return p;
		}
	} else if (vacantSites.size() > 0) { 
		int vacantElemIndex = random.nextInt(vacantSites.size());
		int[] p = (int[])vacantSites.get(vacantElemIndex);
		return (int[])p;
	} else {
		int[] p = new int[2];
		p[0] = rI; p[1] = rJ; 
		return p;
	}
}

public int[] determineLocEpithelial (int cell, int rI, int rJ) {
	LinkedList vacantSites = new LinkedList();
	for (int dI = -1; dI <= 1; dI++) {
		if (((rI + dI) >= 0) || ((rI + dI) <= size - 1)) {
			for (int dJ = -1; dJ <= 1; dJ++) {
				if ((rJ + dJ >= 0) || ((rJ + dJ) <= size - 1)) {
					try {
						if (Cells[rI + dI][rJ + dJ] == 0) 
						{
							if ((EpithelialCells[rI + dI][rJ + dJ]==1) || (EpithelialCells[rI + dI][rJ + dJ]==2)) {
								int[] p = new int[2];
								p[0] = rI + dI; p[1] = rJ + dJ;
								vacantSites.add(p);
							}
						} 
					} catch (Exception ex) {}

				}
			}
		}

	}
	if (vacantSites.size() > 0) {
		int vacantElemIndex = random.nextInt(vacantSites.size());
		int[] p = (int[])vacantSites.get(vacantElemIndex);
		return (int[])p;
	}
	else {
		int[] p = new int[2];
		p[0] = rI; p[1] = rJ; 
		return p;
	}
}

int numberEpithelialEmptySites (int rI, int rJ)
{
	int epithelialSites = 0;
  	for (int dI = -1; dI <= 1; dI++)
   		if (((rI + dI) >= 0) || ((rI + dI) <= size - 1)) 
    			for (int dJ = -1; dJ <= 1; dJ++)
     				if ((rJ + dJ >= 0) || ((rJ + dJ) <= size - 1)) 
						if ((EpithelialCells[rI + dI][rJ + dJ] == 1) || (EpithelialCells[rI + dI][rJ + dJ] == 1)) 
						if (Cells[rI + dI][rJ + dJ]==0) epithelialSites++; 
	return epithelialSites;
}

public boolean isProliferating(int cell, int rI, int rJ) 
    {
        if ((cell==1) && (numberEpithelialEmptySites(rI,rJ)>0) && (random.nextFloat()>TGFB[rI][rJ])) { 
            Age[rI][rJ]=0;
            int[] daughterP = determineLocEpithelial (cell, rI, rJ); 
            Cells[daughterP[0]][daughterP[1]] = EpithelialCells[daughterP[0]][daughterP[1]];
            Age[daughterP[0]][daughterP[1]] = 0.0f;
            TGFBProduction[daughterP[0]][daughterP[1]]=TGFBProduction[rI][rJ];  
            MDEProduction[daughterP[0]][daughterP[1]]=0; 
            MDEProduction[rI][rJ]=0;  
            return true;
        }

        else if (cell != 3 && BM[rI][rJ]>0.01 && (Age[rI][rJ]<1)) return false;
        else {  
            if (((cell ==1) || (cell ==3) || (cell ==4)) && (random.nextFloat() <= (TGFB[rI][rJ]+0.01f)*(proliferatingProbabilityCoefficient[cell]*Age[rI][rJ]))) {  
                Age[rI][rJ]=0;
                int[] daughterP = determineLoc(cell, rI, rJ,false);
                Cells[daughterP[0]][daughterP[1]] = cell;
                Age[daughterP[0]][daughterP[1]] = 0.0f;
                
                TGFBProduction[daughterP[0]][daughterP[1]]=TGFBProduction[rI][rJ];
                if (random.nextFloat()<mutationRate) {  
                    if (random.nextFloat()<0.5f) TGFBProduction[daughterP[0]][daughterP[1]]+=abnormalLuminalCellTGFBProductionRateConstant*wiggleRoom*random.nextFloat();  
                    else  TGFBProduction[daughterP[0]][daughterP[1]]-=abnormalLuminalCellTGFBProductionRateConstant*wiggleRoom*random.nextFloat();
                }  
               
                if (cell ==3) {
                    
                    MDEProduction[daughterP[0]][daughterP[1]]=MDEProduction[rI][rJ];
                    if (random.nextFloat()<mutationRate) {  
                        if (random.nextFloat()<0.5f) MDEProduction[daughterP[0]][daughterP[1]]+=abnormalLuminalCellMDEProductionRateConstant*wiggleRoom*random.nextFloat();  
                        else  MDEProduction[daughterP[0]][daughterP[1]]-=abnormalLuminalCellMDEProductionRateConstant*wiggleRoom*random.nextFloat();
                    }  
                }
                
                if ( ((daughterP[0]==0) || (daughterP[0]==1000) || (daughterP[1]==0) || (daughterP[1]==1000)) && (cell==3) ) {
                    OutDuctTGF=TGFBProduction[daughterP[0]][daughterP[1]]; 
                    OutDuctMDE=MDEProduction[daughterP[0]][daughterP[1]]; 
                    OutDuctCoorX=daughterP[0];
                    OutDuctCoorY=daughterP[1];
                }            
                if ( (outBasLum[daughterP[0]][daughterP[1]] ==-333) && (cell==3) && (OutDuctTGF==-1)) {
                    OutDuctTGF=TGFBProduction[daughterP[0]][daughterP[1]]; 
                    OutDuctMDE=MDEProduction[daughterP[0]][daughterP[1]]; 
                    OutDuctCoorX=daughterP[0];
                    OutDuctCoorY=daughterP[1];
                }
                if (TGFBProduction[daughterP[0]][daughterP[1]]<0) TGFBProduction[daughterP[0]][daughterP[1]]=0;
                if (MDEProduction[daughterP[0]][daughterP[1]]<0) MDEProduction[daughterP[0]][daughterP[1]]=0;
                
                if ((cell ==3) && ((daughterP[0] !=rI) || (daughterP[1] !=rJ))) totalBirths++;
                
                if (TGFBProduction[daughterP[0]][daughterP[1]]>1) TGFBProduction[daughterP[0]][daughterP[1]]=1;
                if (TGFBProduction[daughterP[0]][daughterP[1]]>maxTGFBProd) maxTGFBProd=(TGFBProduction[daughterP[0]][daughterP[1]]);
                else if (TGFBProduction[daughterP[0]][daughterP[1]]<minTGFBProd) minTGFBProd=(TGFBProduction[daughterP[0]][daughterP[1]]);
                newCellTGFBProd+=TGFBProduction[daughterP[0]][daughterP[1]];  
                
                if (MDEProduction[daughterP[0]][daughterP[1]]>1) MDEProduction[daughterP[0]][daughterP[1]]=1;
                if (MDEProduction[daughterP[0]][daughterP[1]]>maxMDEProd) maxMDEProd=(MDEProduction[daughterP[0]][daughterP[1]]);
                else if (MDEProduction[daughterP[0]][daughterP[1]]<minMDEProd) minMDEProd=(MDEProduction[daughterP[0]][daughterP[1]]);
                
                newCellMDEProd+=MDEProduction[daughterP[0]][daughterP[1]];
                
                return true;
            } else return false; 
        }
    }

	
public void iterateTGFB()
{
        float[][] newTGFB = new float[size][size];
        for (int rI = 0; rI < size; rI++)
            for (int rJ = 0; rJ < size; rJ++) {
                int[] top = convertCoordinates(rI, rJ+1);
                int[] left = convertCoordinates(rI-1, rJ);
                int[] right = convertCoordinates(rI+1, rJ);
                int[] south = convertCoordinates(rI, rJ-1);
                
                newTGFB[rI][rJ]
				= TGFB[rI][rJ] +kDt *
				(((TGFB[left[0]][left[1]]
                   + TGFB[south[0]][south[1]]
                   + TGFB[top[0]][top[1]]
                   + TGFB[right[0]][right[1]]
                   - 4.0f * TGFB[rI][rJ])) * (1.0f - BM[rI][rJ]));
                if (BM[rI][rJ]<1.0f) {   
                    newTGFB[rI][rJ] = newTGFB[rI][rJ] - kDt * (1.0f / (4.0f) * 
                    ((BM[right[0]][right[1]] - BM[left[0]][left[1]]) * (TGFB[right[0]][right[1]] - TGFB[left[0]][left[1]])
                    + (BM[top[0]][top[1]] - BM[south[0]][south[1]]) * (TGFB[top[0]][top[1]] - TGFB[south[0]][south[1]])));
                }
                if (Cells[rI][rJ]==6) Cells[rI][rJ]=4;  
                
                switch (Cells[rI][rJ]) {
                    case 1:	
                        if (TGFB[rI][rJ] < 1.0f) 
                            newTGFB[rI][rJ] = newTGFB[rI][rJ]+ basalCellTGFBproductionRate;
						break;
                    case 2:	
                        newTGFB[rI][rJ] =newTGFB[rI][rJ]-luminalCellTGFBconsumptionRateConstant * newTGFB[rI][rJ];
                        break;
                    case 3: 
                        newTGFB[rI][rJ] = newTGFB[rI][rJ] + abnormalLuminalCellTGFBProductionRateConstant*Math.min (1.0f,Age[rI][rJ]/10)*TGFBProduction[rI][rJ];
                        break;
                    case 4: 
                        newTGFB[rI][rJ] = newTGFB[rI][rJ] - stromalTGFBConsumptionRate *TGFB[rI][rJ];
                        if (newTGFB[rI][rJ]>stromalTGFBswitch) {  
                            newTGFB[rI][rJ]=newTGFB[rI][rJ] + stromalCellTGFBproductionRate * TGFB[rI][rJ];
                            Cells[rI][rJ]=6;  
                        }
                        break; 

                }
                newTGFB[rI][rJ] = newTGFB[rI][rJ] * (1-tgfbDecayRate) ; 
                newTGFB[rI][rJ] = newTGFB[rI][rJ] - ECMTGFBconsumption*TGFB[rI][rJ]*BM[rI][rJ]; 
                if (newTGFB[rI][rJ]<0.0f) newTGFB[rI][rJ]=0.0f;
                                    
                if (newTGFB[rI][rJ]>1.0f) newTGFB[rI][rJ]=1.0f;
                
            }
        TGFB = newTGFB;
}

public void iterateMDE()
{
	float[][] newMDE = new float[size][size];
	for (int rI = 0; rI < size; rI++)
		for (int rJ = 0; rJ < size; rJ++) {
			int[] top = convertCoordinates(rI - 1, rJ);
			int[] left = convertCoordinates(rI, rJ - 1);
			int[] right = convertCoordinates(rI, rJ + 1);
			int[] below = convertCoordinates(rI + 1, rJ);
			newMDE[rI][rJ]
				= MDE[rI][rJ] + (kDe *
				(MDE[top[0]][top[1]]
				+ MDE[left[0]][left[1]]
				+ MDE[right[0]][right[1]]
				+ MDE[below[0]][below[1]]
				- 4.0f * MDE[rI][rJ]));

			switch (Cells[rI][rJ]) {  
				case 3:
                    newMDE[rI][rJ]=newMDE[rI][rJ]+abnormalLuminalCellMDEProductionRateConstant*Math.min (1.0f,Age[rI][rJ]/10)*MDEProduction[rI][rJ];
                        break;
			}
            
			newMDE[rI][rJ] = newMDE[rI][rJ]*(1.0f-mdeDecayRate);
            if (newMDE[rI][rJ]<0.0f) newMDE[rI][rJ]=0.0f;
                            
            if (newMDE[rI][rJ]>1.0f) newMDE[rI][rJ]=1.0f;
		}
	MDE = newMDE;
}

public void iterateBM()
{
	float[][] newBM = new float[size][size];
	float kDe = 0.03f;
	for (int rI = 0; rI < size; rI++)
		for (int rJ = 0; rJ < size; rJ++) {
			int[] top = convertCoordinates(rI - 1, rJ);
			int[] left = convertCoordinates(rI, rJ - 1);
			int[] right = convertCoordinates(rI, rJ + 1);
			int[] below = convertCoordinates(rI + 1, rJ);

			newBM[rI][rJ] = BM[rI][rJ] * (1.0f-MDE[rI][rJ]*bmDecayRate);
			switch (Cells[rI][rJ]) {
				case 1:
					newBM[rI][rJ] = newBM[rI][rJ]+basalCellBMProductionRateConstant;
					break;
				case 5: 
					newBM[rI][rJ] = newBM[rI][rJ] + stromalMembraneProductionRate *TGFB[rI][rJ];
					break;
			}
			if (newBM[rI][rJ]>1.0f) newBM[rI][rJ]=1.0f;
			else if (newBM[rI][rJ]<0.0f) newBM[rI][rJ]=0.0f;
		}

	BM = newBM;
}       

	public float [][] getTGFB() { return TGFB; }
	public int [][] getCells () { return Cells; }
	public float [][] getMDE() { return MDE;}
	public float [][] getBM () { return BM;}
	public float [][] getTGFBProd () { return TGFBProduction;}
	public float [][]getMDEProd () { return MDEProduction;}
};
