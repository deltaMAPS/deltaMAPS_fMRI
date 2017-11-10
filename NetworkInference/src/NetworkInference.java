package braininferencesurfacenet;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class NetworkInference {
    
    Vector<Area> areas;
    int dimT, maxLag;
    double significance;
    
    int lCortexNodes;
    int rCortexNodes;
    
    
    public NetworkInference(BrainInferenceSurfaceNet parent)
    {
        this.areas = parent.areas;
        this.dimT = parent.dimT;
        this.maxLag = parent.maxLag;
        this.significance = parent.significanceLevel;
        this.lCortexNodes = parent.lcTimeSeries.length;
        this.rCortexNodes = parent.rcTimeSeries.length;
        exportLargeAreas(250);
   
        networkInference();
        //export weight based link maps
        setAreaStrength(0);
        ComparatorAreaStrength cpStr = new ComparatorAreaStrength();
        java.util.Collections.sort(areas,cpStr);
        printSrengthMap();
        printDistributions();
        ComparatorAreaSize cpSize = new ComparatorAreaSize();
        java.util.Collections.sort(areas,cpSize);
        for(int i = 0; i < 5; i++)
            printLinkMap(areas.get(i), 0);
      
        //exportToInfoMap();
        exportEdgeLists();
        exportAreaInfo();
    }
  
    private void setAreaStrength(int strType)
    {
        for(int i = 0; i < areas.size(); i++)
        {
            double strength = 0;
            Area area = areas.get(i);
            for(int j = 0; j < area.edges.size(); j++)
            {
                if(strType == 0)
                    strength+=Math.abs(area.edges.get(j).weight);
                else
                    strength+=Math.abs(area.edges.get(j).corr);
            }
            areas.get(i).strength = strength;
        }
    }

    private void networkInference()
    {
        System.out.println("Network Inference Started.");
        System.out.println("Significance level: "+this.significance);
        System.out.println("Max lag: "+this.maxLag);
        int nAreas = areas.size();
        int nEdgesTotal = 0, nEdgesDirected = 0;
        for(int i = 0; i < nAreas; i++) //search for an edge between all pairs of areas.
        {
            Area alpha = areas.get(i);
            for(int j = (i+1); j < nAreas; j++)
            {
                Area beta = areas.get(j);            
                
                //step one is to get the correlogram
                double[] correlogram = getCorrelogram(alpha.timeseries, beta.timeseries);
                //step two is to identify significant correlations, use the barlett formula
                double[] significantCorrelations = new double[correlogram.length];
                //(1) get the central part of the time series
                double ts1[] = new double[alpha.timeseries.length-2*maxLag];
                double ts2[] = new double[beta.timeseries.length-2*maxLag];
    
                for(int x = maxLag; x < maxLag+ts1.length; x++)
                {
                    ts1[x-maxLag] = alpha.timeseries[x];
                    ts2[x-maxLag] = beta.timeseries[x];
                }
                //estimate the variance according to Bartlett
                int T = ts1.length;
                double[] autocorrTS1 = new double[2*T+1];
                double[] autocorrTS2 = new double[2*T+1];
                for(int x = 0; x < T; x++)
                {
                    autocorrTS1[T+x] = autocorr(ts1, x);
                    autocorrTS2[T+x] = autocorr(ts2, x);
                    autocorrTS1[T-x] = autocorrTS1[T+x];
                    autocorrTS2[T-x] = autocorrTS2[T+x];
                }
                double var = 0; //after this point var is NOT normalized by T - \tau. 
                for(int x = 0; x < autocorrTS1.length; x++)
                    var+=autocorrTS1[x]*autocorrTS2[x];
                //and now find the significant correlations
                for(int x = 0; x < correlogram.length; x++)
                {
                    int currentLag = Math.abs(maxLag-x);
                    double corr = Math.abs(correlogram[x]); //by taking the absolute corr value is like doing a two-tailed t-test                    
                    //normalize with Bartlett variance
                    corr = corr/Math.sqrt(var/(T-currentLag)); //note: here we further normalize by T-\tau.
                    double pval = 1-DistLib.normal.cumulative(corr, 0.0, 1.0);                
                    if(pval < significance)
                        significantCorrelations[x] = correlogram[x];//passed the test                                        
                }
                AreaEdge e = edgeInference(significantCorrelations);
                if(e!=null)
                {
                    testForLags(significantCorrelations);
                    nEdgesTotal++;
                    //then we have an edge, we need to get the direction and set the weight.
                    double alphaSTD = Math.sqrt(getVariance(ts1, getMean(ts1)));
                    double betaSTD = Math.sqrt(getVariance(ts2, getMean(ts2)));
                    e.weight = alphaSTD*betaSTD*e.corr;
   
                    if(e.lag > 0 ) ///lags are positive alpha -> beta
                    {
                        e.from = alpha.id;
                        e.to = beta.id;
                        areas.get(i).edges.add(e);
                        nEdgesDirected++;
                    }
                    else if(e.lag < 0) //negative lag beta -> alpha
                    {
                        e.from = beta.id;
                        e.to = alpha.id;
                        areas.get(j).edges.add(e);
                        nEdgesDirected++;
                    }
                    else//undirected (lag range includes 0)
                    {
                        AreaEdge e2 = new AreaEdge();
                        e2.corr = e.corr;
                        e2.lag = e.lag;
                        
                        e2.weight = e.weight;
                        e2.undirected=true;
                        e.undirected=true;
                        e.from = alpha.id;
                        e.to = beta.id;
                        areas.get(i).edges.add(e);
                        e2.from = beta.id;
                        e2.to = alpha.id;
                        areas.get(j).edges.add(e2);
                    }
                }
            }
        }
        System.out.println("#Edges (total): "+nEdgesTotal);
        System.out.println("#Edges (directed): "+nEdgesDirected);
    }
    
    private void testForLags(double[] correlogram)
    {
        double maxCorr = 0;
        int lag = 0;
        double zeroLagCorr = correlogram[maxLag];
        
        for(int i = 0; i < correlogram.length; i++)
        {
            if(Math.abs(correlogram[i]) > Math.abs(maxCorr))
            {
                maxCorr = correlogram[i];
                lag = i;
            }
        }
        lag = lag-maxLag;
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("testCorrsNLags.txt",true));
            out.println(maxCorr+"\t"+zeroLagCorr+"\t"+lag);
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    
    
    private AreaEdge edgeInference(double[] significantCorrelations)
    {
        double maxCorr = 0;
        int bestLag = 0;
        //step 1. find the maximum significant correlation in absolute sense
        for(int i = 0; i < significantCorrelations.length; i++)
            if(Math.abs(significantCorrelations[i]) > Math.abs(maxCorr))
            {
                maxCorr = significantCorrelations[i];
                bestLag = i;
            }

        if(maxCorr == 0)//no edge
        {
            return null;
        }
        else
        {
            AreaEdge e = new AreaEdge();
            e.corr = maxCorr;
            e.lag = bestLag-this.maxLag;
            return e;
        }
    }
    
     /*
    returns the autocorrelation of ts at the specified lag
    */
    private double autocorr(double[] tsFoo, int lag)
    {
        double[] ts = Arrays.copyOf(tsFoo, tsFoo.length);
        double mean = getMean(ts);
        double var = getVariance(ts, mean);
        double auto = 0.0;
        for(int i = 0; i < (ts.length-lag); i++)        
            auto+= (ts[i]-mean)*(ts[i+lag]-mean);
        auto = auto/( ((double)ts.length)*var );
        return auto;
    }
    
    
        
    private double[] getCorrelogram(double[] ts1Foo, double[] ts2Foo)
    {                        
        double[] ts1 = Arrays.copyOf(ts1Foo, ts1Foo.length);
        double[] ts2 = Arrays.copyOf(ts2Foo, ts2Foo.length);
        double cijtauArray[] = new double[2*this.maxLag+1];        
     
        //corelation at lag zero
        cijtauArray[this.maxLag]  = pearsonCorrel(ts1, ts2, maxLag,0);        
                
        //start by shifting the second time series to the right
        for(int i = 1; i <= this.maxLag; i++)                
            cijtauArray[this.maxLag+i] = pearsonCorrel(ts1, ts2,maxLag,i);            
        //continue by shifting the time serises to the left
        for(int i = -1; i >= -this.maxLag; i--)                    
            cijtauArray[this.maxLag+i] =  pearsonCorrel(ts1, ts2,maxLag,i);                    
        return cijtauArray;
    }
    
    private double pearsonCorrel(double[] scores1, double[] scores2, int startPos, int lag)
    {    
        double correl = 0.0;
    
        double[] ts1 = new double[scores1.length-2*startPos];
        double[] ts2 = new double[scores1.length-2*startPos];
        
        //get the part of the time series that you want
        for(int i = startPos; i < startPos+ts1.length; i++)
        {
            ts1[i-startPos] = scores1[i];
            ts2[i-startPos] = scores2[i+lag];
        }
        
        //set to unit variance, zero mean first
        ts1 = normalizeToZeroMeanUnitVar(ts1);
        ts2 = normalizeToZeroMeanUnitVar(ts2);
        for(int i = 0; i < ts1.length; i++)
            correl+=ts1[i]*ts2[i];
        
        return correl/(double)ts1.length;        
    }     

    
    private double[] normalizeToZeroMeanUnitVar(double[] ts)
    {
        double mean = getMean(ts);
        double std = Math.sqrt(getVariance(ts, mean));
        for(int i = 0; i < ts.length; i++)
            ts[i] = (ts[i]-mean)/std;
        return ts;
    }
    
    private double getVariance(double[] ts, double mean)
    {
        double var = 0.0;
        for(int i = 0; i < ts.length; i++)        
            var += Math.pow(ts[i]-mean, 2);
        var/=(double)(ts.length-1);
        return var;
    }

    private double getMean(double[] ts)
    {
        double mean = 0.0;
        for(int i = 0; i < ts.length; i++)
            mean+=ts[i];
        mean/=(double)(ts.length);
        return mean;
    }
    
    private double fisherz(double corr)
    {
        return 0.5*Math.log((1+corr)/(1-corr));
    }

    private void printSrengthMap()
    {
        double lCortexStrength[] = new double[lCortexNodes];
        double rCortexStrength[] = new double[rCortexNodes];
        for(int i = areas.size()-1; i >= 0; i --)
        {
            Area area = areas.get(i);
            if(area.cortex == 0)
            {
                for(int j = 0; j < area.surfVoxels.size(); j++)//in case of overlaps show stronger area.
                    if(lCortexStrength[area.surfVoxels.get(j)] < area.strength)
                        lCortexStrength[area.surfVoxels.get(j)] = area.strength;
            }
            else
            {
                for(int j = 0; j < area.surfVoxels.size(); j++)
                    if(rCortexStrength[area.surfVoxels.get(j)] < area.strength)
                        rCortexStrength[area.surfVoxels.get(j)] = area.strength;
            }
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("Lstrength.txt"));
            for(int i = 0; i < lCortexStrength.length; i++)
                out.println(lCortexStrength[i]);
            out.close();
            out = new PrintWriter(new FileWriter("Rstrength.txt"));
            for(int i = 0; i < rCortexStrength.length; i++)
                out.println(rCortexStrength[i]);
            out.close();
            
        }
        catch(IOException ioe)
        {}
    }
    
    private void printDistributions()
    {
        HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
        try
        {
            PrintWriter out1 = new PrintWriter(new FileWriter("strengthDistribution.txt"));
            PrintWriter out2 = new PrintWriter(new FileWriter("stdDistribution.txt"));
            for(int i = 0; i < areas.size(); i++)
            {
                Area area = areas.get(i);
                edges.addAll(area.edges);
                out1.println(area.id+" "+area.strength);
                out2.println(area.id+" "+ Math.sqrt(getVariance(area.timeseries, getMean(area.timeseries))));
            }
            out1.close();
            out2.close();
            out1 = new PrintWriter(new FileWriter("correlationDistribution.txt"));
            out2 = new PrintWriter(new FileWriter("linkWegihtDistribution.txt"));
            Iterator<AreaEdge> it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                out1.println(e.corr);
                out2.println(e.weight);
            }
            out1.close();
            out2.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void printLinkMap(Area area, int linkMapType)
    {
        double[] lCortexLinks = new double[lCortexNodes];
        double[] rCortexLinks = new double[rCortexNodes];
        double[] areaPos = new double[lCortexNodes];
        
        for(int i = 0; i < area.surfVoxels.size(); i++)
            areaPos[area.surfVoxels.get(i)] = -1000000;
        HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
        for(int i = 0; i < areas.size(); i++)
            edges.addAll(areas.get(i).edges);
        Iterator<AreaEdge> it = edges.iterator();
        while(it.hasNext())
        {
            AreaEdge e = it.next();
            Area foo = null;
            if(e.from == area.id)//edge pointing from you to e.to
            {
                foo = new Area(e.to);
                foo = areas.get(areas.indexOf(foo));
            }
            else if(e.to == area.id)//edge pointing from e.from to you
            {
                foo = new Area(e.from);
                foo = areas.get(areas.indexOf(foo));
            }            
            if(foo!= null)
            {
                for(int i = 0; i < foo.surfVoxels.size(); i++)
                {
                    if(foo.cortex == 0)
                    {
                        if(linkMapType == 0)
                            lCortexLinks[foo.surfVoxels.get(i)] = e.weight;
                        else
                            lCortexLinks[foo.surfVoxels.get(i)] = e.corr;
                    }
                    else
                    {
                        if(linkMapType == 0)
                            rCortexLinks[foo.surfVoxels.get(i)] = e.weight;
                        else
                            rCortexLinks[foo.surfVoxels.get(i)] = e.corr;
                    }
                }
            }
            
        }
        
        try
        {
            if(linkMapType == 0)
            {
                PrintWriter out = new PrintWriter(new FileWriter("LlinkMap_"+area.id+".txt"));
                for(int i = 0; i < lCortexLinks.length; i++)
                    out.println(lCortexLinks[i]);
                out.close();
                out = new PrintWriter(new FileWriter("RlinkMap_"+area.id+".txt"));
                for(int i = 0; i < rCortexLinks.length; i++)
                    out.println(rCortexLinks[i]);
                out.close();
                out = new PrintWriter(new FileWriter("linkMapFrom"+area.id+"_"+area.cortex+".txt"));
                for(int i = 0; i < areaPos.length; i++)
                    out.println(areaPos[i]);
                out.close();
            }
            else
            {
                PrintWriter out = new PrintWriter(new FileWriter("LlinkMapCorr_"+area.id+".txt"));
                for(int i = 0; i < lCortexLinks.length; i++)
                    out.println(lCortexLinks[i]);
                out.close();
                out = new PrintWriter(new FileWriter("RlinkMapCorr_"+area.id+".txt"));
                for(int i = 0; i < rCortexLinks.length; i++)
                    out.println(rCortexLinks[i]);
                out.close();
                out = new PrintWriter(new FileWriter("linkMapFromCorr"+area.id+"_"+area.cortex+".txt"));
                for(int i = 0; i < areaPos.length; i++)
                    out.println(areaPos[i]);
                out.close();   
            }
        }
        catch(IOException ioe)
        {}
 
    }
    
    private void exportEdgeLists()
    {
        HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
        for(int i = 0; i < areas.size(); i++)
            edges.addAll(areas.get(i).edges);        
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("edgeList.txt"));
            Iterator<AreaEdge> it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                out.println(e.from+" "+e.to+" "+e.weight);
            }
            out.close();
            out = new PrintWriter(new FileWriter("edgeListPositive.txt"));
            it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                if(e.weight > 0)
                    out.println(e.from+" "+e.to+" "+e.weight);
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void exportAreaInfo()
    {
        try
        {
            PrintWriter outL = new PrintWriter(new FileWriter("areaInfoLeft.txt"));
            PrintWriter outR = new PrintWriter(new FileWriter("areaInfoRight.txt"));
            int lCortexAreas = 0, rCortexAreas = 0;
            for(int i = 0; i < areas.size(); i++)
            {
                Area area = areas.get(i);
                if(area.cortex == 0)
                    lCortexAreas++;
                else
                    rCortexAreas++;
            }
            outL.println("Left cortex, voxels : "+this.lCortexNodes+" areas : "+lCortexAreas);
            outR.println("Right cortex, voxels : "+this.rCortexNodes+" areas : "+rCortexAreas);
            for(int i = 0; i < areas.size(); i++)
            {
                Area area = areas.get(i);
                if(area.cortex == 0)//left cortex.
                {
                    outL.println(area.id);
                    for(int j = 0; j < area.surfVoxels.size(); j++)
                        outL.print(area.surfVoxels.get(j)+" ");
                    outL.print("\r\n");
                }
                else
                {
                    outR.println(area.id);
                    for(int j = 0; j < area.surfVoxels.size(); j++)
                        outR.print(area.surfVoxels.get(j)+" ");
                    outR.print("\r\n");
                }
            }
            outL.close();
            outR.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void exportLargeAreas(int minSize)
    {
        double[] lCortexAreas = new double[this.lCortexNodes];
        double[] rCortexAreas = new double[this.rCortexNodes];
        int largeAreas = 0;
        for(int i = 0; i < areas.size(); i++)
        {
            Area area = areas.get(i);
            if(area.size >= minSize)
            {
                largeAreas++;
                double color = 10*Math.random();
                if(area.cortex == 0)
                {
                    for(int j = 0; j < area.surfVoxels.size(); j++)
                    {
                        int voxelID = area.surfVoxels.get(j);
                        double areaColor = color;
                        if(lCortexAreas[voxelID]!= 0)//overlap, make color negative
                            areaColor = -2.0;
                        else
                            areaColor = color;
                        lCortexAreas[voxelID] = areaColor;
                    }
                }
                else
                {
                    for(int j = 0; j < area.surfVoxels.size(); j++)
                    {
                        int voxelID = area.surfVoxels.get(j);
                        double areaColor = color;
                        if(rCortexAreas[voxelID]!= 0)//overlap, make color negative
                            areaColor = -2;
                        else
                            areaColor = color;
                        rCortexAreas[voxelID] = areaColor;
                    }
                }
            }
        }
        try
        {
            PrintWriter outL = new PrintWriter(new FileWriter("largeAreasL.txt"));
            PrintWriter outR = new PrintWriter(new FileWriter("largeAreasR.txt"));
            for(int i = 0; i < rCortexAreas.length; i++)
                outR.println(rCortexAreas[i]);
            for(int i = 0; i < lCortexAreas.length; i++)
                outL.println(lCortexAreas[i]);
            outL.close();
            outR.close();
        }
        catch(IOException ioe)
        {}
        System.out.println("Exported "+largeAreas+" large areas.");
    }
    
    private void exportToInfoMap()
    {
        //Step 1. Get all unique edges in a HashSet
        HashSet<AreaEdge> edges = new HashSet<AreaEdge>();
        for(int i = 0; i < areas.size(); i++)
            edges.addAll(areas.get(i).edges);
        //Step 2. Map all area IDs to [1,N]
        int newID = 1;
        Hashtable<Integer,Integer> mapAreaIDto1N = new Hashtable<Integer,Integer>();
        Iterator<AreaEdge> it = edges.iterator();
        while(it.hasNext())
        {
            AreaEdge e = it.next();
            if(e.weight > 0)//infomap can not handle negative weights.
            {
                if(!mapAreaIDto1N.containsKey(e.from))
                {
                    mapAreaIDto1N.put(e.from, newID);
                    newID++;
                }
                if(!mapAreaIDto1N.containsKey(e.to))
                {
                    mapAreaIDto1N.put(e.to, newID);
                    newID++;                        
                }
            }
        }
        
        try
        {
            //Step 3. Export edge list
            PrintWriter out = new PrintWriter(new FileWriter("infomapEdgeList.txt"));
            it = edges.iterator();
            while(it.hasNext())
            {
                AreaEdge e = it.next();
                if(e.weight > 0)
                {
                    int newIDfrom = mapAreaIDto1N.get(e.from);
                    int newIDto = mapAreaIDto1N.get(e.to);
                    out.println(newIDfrom+" "+newIDto+" "+e.weight);
                }
            }
            out.close();
            //Step 4. Export map of area ids
            out = new PrintWriter(new FileWriter("infomapIDs.txt"));
            java.util.Enumeration<Integer> oldIDs = mapAreaIDto1N.keys();
            java.util.Enumeration<Integer> infoIDs = mapAreaIDto1N.elements();
            while(oldIDs.hasMoreElements())            
                out.println(oldIDs.nextElement()+" "+infoIDs.nextElement());            
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
}

class ComparatorAreaStrength implements Comparator<Area>
{
    public int compare(Area alpha, Area beta)
    {
        if(alpha.strength> beta.strength)
            return -1;
        else if(alpha.strength < beta.strength)
            return 1;
        else return 0;
    }
}

class ComparatorAreaSize implements Comparator<Area>
{
    public int compare(Area alpha, Area beta)
    {
        if(alpha.surfVoxels.size() > beta.surfVoxels.size())
            return -1;
        else if(alpha.surfVoxels.size() < beta.surfVoxels.size())
            return 1;
        else return 0;
    }
}
