package braininferencesurface;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class NetworkInference {
    
    boolean reportExactLag = false;
    int dimT;
    int maxLag;
    double netSigLevel;
    Vector<Area> areas;
    double[][] timeseries;
    
    
    public NetworkInference(Vector<Area> areas, String timeseriesFile, double netSigLevel, int maxLag)
    {
        this.areas = areas;
        this.netSigLevel = netSigLevel;
        importTimeseries(timeseriesFile);
        this.maxLag = maxLag;
        System.out.println("Constructing area time series.");
        constructAreaTimeSeries();
    }
    
    private void constructNetwork()
    {        
        
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);            
            for(int j = (i+1); j < areas.size(); j++)
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
                //estimate the variance according to Barlett
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
                double var = 0;
                for(int x = 0; x < autocorrTS1.length; x++)
                    var+=autocorrTS1[x]*autocorrTS2[x];
                //and now find the significant correlations
                for(int x = 0; x < correlogram.length; x++)
                {
                    int currentLag = Math.abs(maxLag-x);
                    double corr = Math.abs(correlogram[x]); //by taking the absolute corr value is like doing a two-tailed t-test
                    //fisher transform
                    corr = fisherZ(corr);
                    //normalize
                    corr = corr/Math.sqrt(var/(T-currentLag));
                    double pval = 1-DistLib.normal.cumulative(corr, 0.0, 1.0);                
                    if(pval < netSigLevel)
                        significantCorrelations[x] = correlogram[x];//passed the test                                        
                }                                  
                //ILIA : the code above works.
                //now find the max. significant correlation and the lag range.
                double linkCorr = 0.0;
                int linkLag = 0, minLagRange = -100000, maxLagRange = -100000;
                for(int x = 0; x < significantCorrelations.length; x++)
                {
                    if(Math.abs(significantCorrelations[x]) > Math.abs(linkCorr))
                    {
                        linkCorr = significantCorrelations[x];
                        linkLag = x-maxLag;
                    }
                    if(Math.abs(significantCorrelations[x]) > 0 && minLagRange == -100000)
                        minLagRange = x-maxLag;
                    if(Math.abs(significantCorrelations[x]) > 0)
                        maxLagRange = x-maxLag;
                }
                
                if(linkCorr!= 0)//we found an edge!
                {
                    if(reportExactLag)
                    {    
                        AreaEdge e = null;
                        if(linkLag >= 0)                    
                            e = new AreaEdge(alpha.id, beta.id, 0, linkLag);
                        else
                            e = new AreaEdge(beta.id, alpha.id, 0, linkLag);
                        e.sij=linkCorr;
                        
                        
                        e.minLag = minLagRange;
                        e.maxLag = maxLagRange;
                        double alphaSTD = getVariance(ts1, getMean(ts1));
                        double betaSTD = getVariance(ts2, getMean(ts2));
                        e.weight=Math.sqrt(alphaSTD)*Math.sqrt(betaSTD)*linkCorr;
                        areas.get(i).edges.add(e);
                        areas.get(j).edges.add(e);
                    }
                    if(!reportExactLag)
                    {
                        AreaEdge e = null;
                        //if min lag positive then A-> B
                        if(minLagRange > 0)
                            e = new AreaEdge(alpha.id, beta.id, 0, linkLag);
                        else if(maxLagRange < 0)//if max lag negative B -> A
                            e = new AreaEdge(beta.id, alpha.id, 0, linkLag);
                        else if(minLagRange <= 0 && maxLagRange >= 0) //then bi-directed link
                            e = new AreaEdge(beta.id, alpha.id, 0, 0);
                        try
                        {
                            e.sij=linkCorr;
                        }
                        catch(NullPointerException npe)
                        {
                            System.out.println();
                        }
                        e.minLag = minLagRange;
                        e.maxLag = maxLagRange;
                        double alphaSTD = getVariance(ts1, getMean(ts1));
                        double betaSTD = getVariance(ts2, getMean(ts2));
                        e.weight=Math.sqrt(alphaSTD)*Math.sqrt(betaSTD)*linkCorr;
                        areas.get(i).edges.add(e);
                        areas.get(j).edges.add(e);
                    }
                }                                                                                
            }
        }
        
        int edgesin = 0, edgesout = 0;
        for(int i = 0; i < areas.size(); i++)
        {
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.edges.size(); j++)
            {
                AreaEdge e = alpha.edges.get(j);                
                if(e.from == alpha.id)
                    edgesout++;
                else edgesin++;
            }
        }
        System.out.println("Edges: "+(edgesin+edgesout));
        System.out.println("Edges in: "+edgesin);
        System.out.println("Edges out: "+edgesout);
        double density = (edgesin+edgesout)/(double)((areas.size()*areas.size()-1)/2);
        System.out.println("Density: "+density);
        
    }

    private double fisherZ(double corr)
    {
        if(corr == 1)
            corr = 0.999999999;
         return 0.5*Math.log( (1+corr)/(1-corr));
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


    private void constructAreaTimeSeries()
    {
        int nAreas = areas.size();
        for(int i = 0; i < nAreas; i++)
        {
            double nNodes = 0;
            double areaTS[] = new double[dimT];
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.nodes.size(); j++)//add the time series of all grid cells.
            {
                int nodeid = alpha.nodes.get(j);
                nNodes++;
                double[] ts = timeseries[nodeid];
                for(int t = 0; t < dimT; t++)
                    areaTS[t]+=ts[t];
            }
            //and take the average
            for(int t = 0; t < dimT; t++)
                areaTS[t]=areaTS[t]/nNodes;
            areas.get(i).timeseries=areaTS;
        }
    }
    
    //DEBUG STATUS: FALSE
    private void importTimeseries(String timeseriesFile)
    {
        //first of all the time series are in a file in which each line corresponds to the ith nodes time series
        //so before adding anything we need to make sure that the nodes are sorted by id in ascending order 
        System.out.println("Loading time series.");
        try
        {
            int nLines = 0;
            BufferedReader in = new BufferedReader(new FileReader(timeseriesFile));
            String line = in.readLine();
            String[] args = line.split(",");
            this.dimT = args.length;
            try
            {
                while(line!= null)
                {
                    nLines++;
                    line = in.readLine();                    
                }
            }
            catch(NullPointerException npe)
            {}
            timeseries = new double[nLines][dimT];
            in = new BufferedReader(new FileReader(timeseriesFile));
            line = in.readLine();
            int cc = 0;
            try
            {
                while(line!= null)
                {
                    args = line.split(",");
                    for(int i = 0; i < dimT; i++)
                        timeseries[cc][i] = Double.parseDouble(args[i]);
                    cc++;
                    line = in.readLine();
                }
            }
            catch(NullPointerException npe)
            {}
            in.close();
        }
        catch(IOException ioe)
        {}
        System.out.println("Loaded "+timeseries.length+" time series. dimT: "+timeseries[0].length);
    }
    
    
}
