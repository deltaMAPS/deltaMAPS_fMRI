package braininferencesurface;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class PeakIdentification {
    Vector<SurfNode> surfNodes;
    double areaThreshold;
    int neighSize, dimT;

    public PeakIdentification(Vector<SurfNode> surfNodes, int neighSize, double areaThreshold) {
        System.out.println("Peak identification started.");
        this.surfNodes = surfNodes;
        this.dimT = surfNodes.firstElement().timeseries.length;
        this.areaThreshold = areaThreshold;
        this.neighSize = neighSize;
        System.out.println("Neighborhood size set to: "+this.neighSize);
        
        constructLocalNeighborhoods();
        inferScores();
        inferPeaks();
        printScores();
        printPeaks();
    }
    
    private void inferPeaks()
    {
        int nPeaks = 0;
        System.out.println("Peak inference started.");
        for(int i = 0; i < surfNodes.size(); i++)
        {
            SurfNode node = surfNodes.get(i);
            if(node.hasData)
            {
                if(node.localNeighScore > this.areaThreshold)
                {
                    boolean isPeak = true;
                    for(int j = 0; j < node.localNeighborhood.size(); j++)
                    {
                        SurfNode node2 = new SurfNode(node.localNeighborhood.get(j));
                        node2 = surfNodes.get(surfNodes.indexOf(node2));
                        if(node2.localNeighborhood == null)
                            isPeak = false;
                        if(node2.localNeighScore > node.localNeighScore)
                            isPeak = false;
                    }
                    surfNodes.get(i).isPeak = isPeak;
                    if(isPeak)
                        nPeaks++;
                }
            }
        }
        System.out.println("Found a total of "+nPeaks+" peaks.");
    }
    
    private void inferScores()
    {
        System.out.println("Score inference started.");
        for(int i = 0; i < surfNodes.size(); i++)
        {
            SurfNode n = surfNodes.get(i);
            if(n.localNeighborhood!=null)
            {
                double score = getScore(n.localNeighborhood);
                surfNodes.get(i).localNeighScore=score;
            }
        }
    }
    
    private double getScore(Vector<Integer> nodes)
    {
        int N = nodes.size();
        double denom = N*(N-1)/2;
        double score = 0.0;
        for(int i = 0; i < N; i++)
        {
            SurfNode n1 = new SurfNode(nodes.get(i));
            n1 = surfNodes.get(surfNodes.indexOf(n1));
            for(int j = (i+1); j < N; j++)
            {
                SurfNode n2 = new SurfNode(nodes.get(j));
                n2 = surfNodes.get(surfNodes.indexOf(n2));
                score+=pearsonCorrel(n1.timeseries, n2.timeseries);
            }
        }
        return score/denom;
    }
    
    private void constructLocalNeighborhoods()
    {
        int nodesWithNeighborhood = 0;
        System.out.println("Constructing local neighborhoods.");
        for(int i = 0; i < surfNodes.size(); i++)
        {
            SurfNode n = surfNodes.get(i);
            if(n.hasData)
            {   
                Vector<Integer> localNeighborhood = getLocalNeighborhood(n.nodeID);
                if(localNeighborhood.size() == neighSize+1)//plus 1 because local neighborhood contains node i
                {
                    surfNodes.get(i).addLocalNeighbors(localNeighborhood);
                    nodesWithNeighborhood++;
                }
                else
                    surfNodes.get(i).localNeighborhood = null;                            
            }
            else surfNodes.get(i).localNeighborhood = null;
        }
        System.out.println("Found "+nodesWithNeighborhood+" nodes with local neighborhood");
    }
    
    private Vector<Integer> getLocalNeighborhood(int nodeid)
    {
        Vector<Integer> localNeighborhood = new Vector<Integer>();
        localNeighborhood.add(nodeid);
        SurfNode startNode = new SurfNode(nodeid);
        startNode = surfNodes.get(surfNodes.indexOf(startNode));
        
        while(localNeighborhood.size() != neighSize+1)
        {
            HashSet<Integer> candidates = new HashSet<Integer>(); //use a hashset to prevent duplicates.
            //get the geographically adjacent nodes to the nodes currently in the local neighborhood.
            for(int i = 0; i < localNeighborhood.size(); i++)
            {
                SurfNode n = new SurfNode(localNeighborhood.get(i));
                n = surfNodes.get(surfNodes.indexOf(n));
                
                
                for(int j = 0; j < n.edges.size(); j++)
                {
                    int neighID = n.edges.get(j);
                    if(!localNeighborhood.contains(neighID)) //a candidate should not be already in the local neighborhood
                        candidates.add(neighID);
                }
            }
            //next if we have zero candidates return
            if(candidates.size() == 0)
                return localNeighborhood;
            //next, find and add the node in the candidates with the maximal correlation to the start point
            double bestCorr = -1;
            int bestNeighID = -1;
            java.util.Iterator<Integer> it = candidates.iterator();
            while(it.hasNext())
            {
                int candidateID = it.next();
                SurfNode cnode = new SurfNode(candidateID);
                cnode = surfNodes.get(surfNodes.indexOf(cnode));
                double corr = pearsonCorrel(startNode.timeseries, cnode.timeseries);
                if(corr > bestCorr)
                {
                    bestCorr = corr;
                    bestNeighID = candidateID;
                }               
            }
            //add him to the local neighborhood
            localNeighborhood.add(bestNeighID);
        }
            
        return localNeighborhood;
    }
    
 
    
    private double pearsonCorrel(double[] ts1, double[] ts2)
    {
       double[] scores1 = Arrays.copyOf(ts1, dimT);
       double[] scores2 = Arrays.copyOf(ts2, dimT);
       scores1 = normalizeToZeroMeanUnitVar(scores1);
       scores2 = normalizeToZeroMeanUnitVar(scores2);
       double correl = 0.0;
       for(int i = 0; i < scores1.length; i++)
           correl+=scores1[i]*scores2[i];
       correl = correl/(double)dimT;
       
       return correl;
    }
    
    private double[] normalizeToZeroMeanUnitVar(double[] ts)
    {
        double mean = getMean(ts);
        double std = Math.sqrt(getVariance(ts, mean));
        for(int i = 0; i < ts.length; i++)
            ts[i] = (ts[i]-mean)/std;
        return ts;
    }
    
    private double getMean(double[] ts)
    {
        double mean = 0.0;
        for(int i = 0; i < ts.length; i++)
            mean+=ts[i];
        mean/=(double)(ts.length);
        return mean;
    }
       
    
    private double getVariance(double[] ts, double mean)
    {
        double var = 0.0;
        for(int i = 0; i < ts.length; i++)        
            var += Math.pow(ts[i]-mean, 2);
        var/=(double)(ts.length-1);
        return var;
    }
    
    private void printScores()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("peakScores.txt"));
            for(int i = 0; i < surfNodes.size(); i++)
            {
                SurfNode node = surfNodes.get(i);
                if(node.localNeighScore == -1)
                    out.println(0);
                else out.println(node.localNeighScore);
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void printPeaks()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("peakPos.txt"));
            for(int i = 0; i < surfNodes.size(); i++)
            {
                SurfNode n = surfNodes.get(i);
                if(n.isPeak)
                    out.println(1);
                else
                    out.println(0);
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    
}
