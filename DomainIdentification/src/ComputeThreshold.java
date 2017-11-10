package braininferencesurface;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class ComputeThreshold {
    
    Vector<SurfNode> surfNodes;
    double significanceLevel;
    double threshold;
    int dimT;
    int sampleSize = 10000;
    
    public ComputeThreshold(Vector<SurfNode> surfNodes, double significanceLevel)
    {
        System.out.println("Threshold inference starts using "+sampleSize+" samples");
        System.out.println("Alpha_A: "+significanceLevel);
        
        this.surfNodes = surfNodes;
        this.dimT = surfNodes.firstElement().timeseries.length;
        this.significanceLevel = significanceLevel;
        exportCorrelationSamples(10000);
        this.threshold = getThresholdBarlett();
        System.out.println("Threshold set to: "+this.threshold);
    }
    
    private double getThresholdBarlett()
    {
        int currentSamples = 0;
        double avgThreshold = 0.0, nCorrs = 0.0;
        double minSignificantCorrel = 1.0;
                    
        while(currentSamples < sampleSize)
        {
            //get two time series
            int pos1 = (int)Math.floor(Math.random()*(surfNodes.size()-1));
            int pos2 = (int)Math.floor(Math.random()*(surfNodes.size()-1));
            
            if(currentSamples%1000 == 0)
                System.out.println("Progress: "+((double)currentSamples/(double)sampleSize)+"%");
            //get their cross correlation
            if(surfNodes.get(pos1).hasData && surfNodes.get(pos2).hasData)
            {   
                double[] ts1 = surfNodes.get(pos1).timeseries;
                double[] ts2 = surfNodes.get(pos2).timeseries;
                
                double corr = pearsonCorrel(ts1, ts2);
                //estimate variance using Barlett
                double var = varianceBarlett(ts1, ts2);
                //fisher transform the correlation
                double corrF = fisherZ(corr);
                //normalize by std of Barlett
                double zijF = corrF/Math.sqrt(var);
                //test significance
                double pval = 1-DistLib.normal.cumulative(zijF, 0.0, 1.0);                
                if(pval < significanceLevel)
                {  //passed
                    avgThreshold+=corr;
                    nCorrs++;        
                    if(minSignificantCorrel > corr)
                        minSignificantCorrel = corr;
                }

                currentSamples++;
            }
        }

        System.out.println("Min significant correlation: "+minSignificantCorrel);
        return avgThreshold/nCorrs;
    }
    
    
    private double varianceBarlett(double[] ts11, double[] ts21)
    {
        double[] ts1 = Arrays.copyOf(ts11, dimT);
        double[] ts2 = Arrays.copyOf(ts21, dimT);
        //autocorrelation from [-T,T] for ts1 and ts2
        double[] autocorrTS1 = new double[2*dimT+1];
        double[] autocorrTS2 = new double[2*dimT+1];
        for(int i = 0; i < dimT; i++)
        {
            autocorrTS1[dimT+i] = autocorr(ts1, i);
            autocorrTS2[dimT+i] = autocorr(ts2, i);
            autocorrTS1[dimT-i] = autocorrTS1[dimT+i];
            autocorrTS2[dimT-i] = autocorrTS2[dimT+i];
        }
        double var = 0.0;
        for(int i = 0; i < autocorrTS1.length; i++)
            var+= autocorrTS1[i]*autocorrTS2[i];
        return var/(double)dimT;        
    }
    
    /*
    returns the autocorrelation of ts at the specified lag
    */
    private double autocorr(double[] ts, int lag)
    {
        double mean = getMean(ts);
        double var = getVariance(ts, mean);
        double auto = 0.0;
        for(int i = 0; i < (dimT-lag); i++)        
            auto+= (ts[i]-mean)*(ts[i+lag]-mean);
        auto = auto/( ((double)dimT)*var );
        return auto;
    }
    
    private double fisherZ(double corr)
    {
        if(corr == 1)
            corr = 0.999999999;
         return 0.5*Math.log( (1+corr)/(1-corr));
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
    
    private void exportCorrelationSamples(int nSamples)
    {
        int nNodes = surfNodes.size();
        try
        {
            //correlations between random nodes        
            PrintWriter out = new PrintWriter(new FileWriter("randCorrSamples.txt"));
            int cc = 0;        
            while(cc != nSamples)
            {
                //select two random surfNodes.
                int pos1 = (int)Math.floor(Math.random()*(nNodes-1));
                int pos2 = (int)Math.floor(Math.random()*(nNodes-1));
                SurfNode sf1 = surfNodes.get(pos1);
                SurfNode sf2 = surfNodes.get(pos2);
                if(sf1.hasData && sf2.hasData)
                {
                    cc++;
                    out.println(pearsonCorrel(sf1.timeseries, sf2.timeseries));
                }
            }
            out.close();
            //correlations between adjacent nodes here
            out = new PrintWriter(new FileWriter("neighCorrSamples.txt"));
            cc = 0;        
            while(cc != nSamples)
            {
                //select two random surfNodes.
                int pos1 = (int)Math.floor(Math.random()*(nNodes-1));
                SurfNode sf1 = surfNodes.get(pos1);
                if(sf1.hasData)
                {
                    cc++;
                    SurfNode sf2 = new SurfNode(sf1.edges.firstElement());
                    sf2 = surfNodes.get(surfNodes.indexOf(sf2));
                    out.println(pearsonCorrel(sf1.timeseries, sf2.timeseries));
                }
            }
            out.close();
        }
        catch(IOException ioe)
        {}    
            
    }

    
}
