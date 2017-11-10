package braininferencesurfacenet;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class BrainInferenceSurfaceNet {

    Vector<Area> areas;
    double[][] lcTimeSeries;
    double[][] rcTimeSeries;
    int areaID = 0;
    
    double significanceLevel;
    int maxLag, dimT;
    
    public static void main(String[] args) {
        System.out.println("Input:");
        System.out.println("1. Areas for left cortex");
        System.out.println("2. Areas for right cortex");
        System.out.println("3. Left Cortex time series file");
        System.out.println("4. Right cortex time series file");
        System.out.println("5. Max lag");
        System.out.println("6. Network significance level.");
        String leftCAreasFile = "LfinalAreas122620.txt";
        String rightCAreasFile = "RfinalAreas122620.txt";
        String leftCTsFile = "Ltimeseries122620.txt";
        String rightCTsFile = "Rtimeseries122620.txt";
        int maxLag = 1;//Integer.parseInt(args[4]);
        double significanceLevel = 0.000001;//Double.parseDouble(args[5]);
        BrainInferenceSurfaceNet bisn =new BrainInferenceSurfaceNet(leftCAreasFile, rightCAreasFile, leftCTsFile, 
                rightCTsFile, maxLag, significanceLevel);
    }
    
     public BrainInferenceSurfaceNet(String leftCAreasFile, String rightCAreasFile, String leftCTsFile,
         String rightCTsFile, int maxLag, double significance)
     {
        this.maxLag = maxLag;
        this.significanceLevel = significance;
        //Step 1. load the areas.
        System.out.println("Loading areas.");
        areas = new Vector<Area>();
        loadAreas(leftCAreasFile, 0);
        loadAreas(rightCAreasFile, 1);
        System.out.println("Loaded "+areas.size()+" areas.");
        
        //Step 2. load the surface voxel timeseries
        System.out.println("Loading time series.");
        loadTimeSeries(leftCTsFile, 0);
        loadTimeSeries(rightCTsFile, 1);
        System.out.println("dimT: "+this.dimT);
        System.out.println("Left cortex time series: "+this.lcTimeSeries.length);
        System.out.println("Right cortex time series: "+this.rcTimeSeries.length);
        //Step 3. construct the area time series
        System.out.println("Constructing area time series.");
        constructAreaTimeSeries();
        System.out.println("Starting Network Inference.");
        NetworkInference ni  = new NetworkInference(this);
        
    }
    
    private void constructAreaTimeSeries()
    {
        for(int i = 0; i < areas.size(); i++)
        {
            Area area = areas.get(i);
            Vector<Integer> surfVoxels = area.surfVoxels;
            double[] areaTs = new double[dimT];
            for(int j = 0; j < surfVoxels.size(); j++)
            {
                if(area.cortex == 0)
                {
                    int voxelID = surfVoxels.get(j);
                    double[] ts = lcTimeSeries[voxelID];
                    for(int t = 0; t < dimT; t++)
                        areaTs[t]+=ts[t];
                }
                else
                {
                    int voxelID = surfVoxels.get(j);
                    double[] ts = rcTimeSeries[voxelID];
                    for(int t = 0; t < dimT; t++)
                        areaTs[t]+=ts[t];
                }
            }
            for(int j = 0; j < dimT; j++)//get the average
                areaTs[j] = areaTs[j]/(double)area.size;
            areas.get(i).addTimeSeries(areaTs);
        }
    }
    
    private void loadTimeSeries(String file, int cortex)
    {
        //first pass in to read the file, get the time dimension and the number
        //of voxels of the left/right cortex
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = in.readLine();
            this.dimT = line.split(",").length;
            int nVoxels = 0;
            try
            {
                while(line!= null)
                {
                    nVoxels++;
                    line = in.readLine();
                }
            }
            catch(IOException ioe)
            {}
            in.close();
            //initialize time series array
            if(cortex == 0)
                this.lcTimeSeries = new double[nVoxels][dimT];
            else
                this.rcTimeSeries = new double[nVoxels][dimT];
            //now add the time series
            in = new BufferedReader(new FileReader(file));
            line = in.readLine();
            try
            {
                int cc = 0;//position to add the time series in the 2d array
                while(line!=null)
                {
                    String args[] = line.split(",");
                    double[] ts = new double[dimT];
                    for(int i = 0; i < args.length; i++)
                        ts[i] = Double.parseDouble(args[i]);
                    if(cortex == 0)
                        lcTimeSeries[cc]=ts;
                    else
                        rcTimeSeries[cc]=ts;
                    cc++;
                    line = in.readLine();
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
        }
        catch(IOException ioe)
        {
            System.out.println("Failed to load: "+file);
            ioe.printStackTrace();
            System.exit(1);
        }
    }
    
    /*
    cortex = 0, left cortex; cortex = 1, right cortex
    */
    private void loadAreas(String file, int cortex)
    {
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = in.readLine();
            try
            {
                while(line!= null)
                {
                    //first line is dummy ignore
                    line = in.readLine();
                    String[] args = line.split("\t");
                    Area area = new Area(areaID);
                    areaID++;
                    area.cortex = cortex;
                    for(int i = 0; i < args.length; i++)
                        area.surfVoxels.add(Integer.parseInt(args[i]));
                    area.size = area.surfVoxels.size();
                    areas.add(area);
                    line = in.readLine();
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
        }
        catch(IOException ioe)
        {
            System.out.println("Failed to load: "+file);
            ioe.printStackTrace();
            System.exit(0);
        }
    }
    
}
