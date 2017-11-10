package braininferencesurface;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class BrainInferenceSurface {

    Vector<SurfNode> surfNodes; //here we will store all surface nodes, their time series and their spatially adjacent nodes.
    double areThreshold;
    int maxLag;
    
    
    public static void main(String[] args) {
        
        System.out.println("INPUT:");
        System.out.println("1. Time series file. Successive values hsould be comma separated.");
        System.out.println("2. Surface face file");
        System.out.println("3. Max lag");
        System.out.println("3. Area Identification significance level");
        System.out.println("4. Neighborhood expansion factor K");
                
        String timeseriesFile = args[0];//"timeseriesOutFuncNoMask.txt";// args[0];
        String faceFile = args[1];//"facesConte.txt";//args[1];
        int maxLag = Integer.parseInt(args[2]);
        double areaSigLevel =  Double.parseDouble(args[3]);
        int neighFactor =  Integer.parseInt(args[4]);
        
        BrainInferenceSurface bis = new BrainInferenceSurface(timeseriesFile, faceFile, maxLag, areaSigLevel, neighFactor);
    }
    
    public BrainInferenceSurface(String timeseriesFile, String faceFile, int maxLag, double areaSigLevel
                                    ,int neighFactor)
    {
        long startTime = System.currentTimeMillis();
        this.maxLag = maxLag;
        importFaceFile(faceFile);
        importTimeseries(timeseriesFile);
        printSpatialDegree(this.surfNodes);
        cleanUpData();
        System.out.println("--------------------Stage 2--------------------");
        ComputeThreshold ct = new ComputeThreshold(surfNodes, areaSigLevel);
        this.areThreshold = ct.threshold;
        System.out.println("--------------------Stage 3--------------------");
        PeakIdentification pf = new PeakIdentification(surfNodes, neighFactor, areThreshold);
        System.out.println("--------------------Stage 4--------------------");
        AreaIdentification ai = new AreaIdentification(pf.surfNodes, areThreshold);
        long endTime = System.currentTimeMillis();
        long duration = (endTime-startTime)*1000;
        System.out.println("Finished in: "+duration+" seconds");
    }
    
    //DEBUG STATUS: TRUE.
    private void cleanUpData()
    {
        System.out.println("Cleaning up data.");
        //some of the nodes do not have any signal (their time series contain zero values)
        //remove these nodes from the analysis by removing them from the "edges" field.
        Vector<Integer> nodeidsWithNoData = new Vector<Integer>();
        for(int i = 0; i < surfNodes.size(); i++)
        {
            SurfNode n = surfNodes.get(i);
            boolean hasData = n.hasData();
            if(!hasData)
            {
                surfNodes.get(i).hasData = false;
                nodeidsWithNoData.add(n.nodeID);
            }
        }
        
       
        for(int i = 0; i < surfNodes.size(); i++)
            for(int j = 0; j < nodeidsWithNoData.size(); j++)
                surfNodes.get(i).removeEdge(nodeidsWithNoData.get(j));
        System.out.println("Found "+nodeidsWithNoData.size()+" nodes with no data.");       
    }

    //DEBUG STATUS: TRUE
    private void importTimeseries(String timeseriesFile)
    {
        //first of all the time series are in a file in which each line corresponds to the ith nodes time series
        //so before adding anything we need to make sure that the nodes are sorted by id in ascending order 
        System.out.println("Adding time series to nodes.");
        NodeIDComparator mycomp = new NodeIDComparator();
        java.util.Collections.sort(surfNodes,mycomp);
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(timeseriesFile));
            String line = in.readLine();
            int currentNode = 0;
            try
            {
                while(line!= null)
                {
                    String args[] = line.split(",");
                    int dimT = args.length;
                    surfNodes.get(currentNode).addAreaTimeseriesFromFile(args, dimT,maxLag);
                    line = in.readLine();
                    currentNode++;
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
        }
        catch(IOException ioe)
        {}
    }
    
    //DEBUG STATUS: TRUE.
    private void importFaceFile(String faceFile)
    {
        surfNodes = new Vector<SurfNode>();
        //first pass, get all unique node ids.
        HashSet<Integer> uniqueIDs = new HashSet<Integer>();
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(faceFile));
            String line = in.readLine();
            try
            {
                while(line!=null)
                {
                    String args[] = line.split(",");
                    for(int j = 0; j < args.length; j++)
                        uniqueIDs.add(Integer.parseInt(args[j])-1);
                    line = in.readLine();
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
            System.out.println("Found "+uniqueIDs.size()+" nodes");
        }
        catch(IOException ioe)
        {}
        System.out.println("Initializing nodes");
        java.util.Iterator<Integer> it = uniqueIDs.iterator();
        while(it.hasNext())
        {
            Integer nodeid = it.next();
            SurfNode n = new SurfNode(nodeid);
            surfNodes.add(n);
        }
        System.out.println("Initializing node adjacencies");
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(faceFile));
            String line = in.readLine();
            try
            {
                while(line!= null)
                {
                    String triplet[] = line.split(",");
                    SurfNode n1 = new SurfNode(Integer.parseInt(triplet[0])-1);
                    SurfNode n2 = new SurfNode(Integer.parseInt(triplet[1])-1);
                    SurfNode n3 = new SurfNode(Integer.parseInt(triplet[2])-1);
                    int indexN1 = surfNodes.indexOf(n1);
                    int indexN2 = surfNodes.indexOf(n2);
                    int indexN3 = surfNodes.indexOf(n3);
                    surfNodes.get(indexN1).addEdge(n2.nodeID);
                    surfNodes.get(indexN1).addEdge(n3.nodeID);
                    surfNodes.get(indexN2).addEdge(n1.nodeID);
                    surfNodes.get(indexN2).addEdge(n3.nodeID);
                    surfNodes.get(indexN3).addEdge(n1.nodeID);
                    surfNodes.get(indexN3).addEdge(n2.nodeID);
  
                    
                    line = in.readLine();
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
        }
        catch(IOException ioe)
        {}
                
    }
    
    private void printSpatialDegree(Vector<SurfNode> surfNodes)    
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("spaptDegree.txt"));
            for(int i = 0; i < surfNodes.size(); i++)
                out.println(surfNodes.get(i).edges.size());
            out.close();
        }
        catch(IOException ioe)
        {}
    }
}


class NodeIDComparator implements Comparator<SurfNode>
{
    @Override
    public int compare(SurfNode n1, SurfNode n2)
    {
        if(n1.nodeID > n2.nodeID)
            return 1;
        else if(n1.nodeID < n2.nodeID)
            return -1;
        else return 0;
                                                
    }
}