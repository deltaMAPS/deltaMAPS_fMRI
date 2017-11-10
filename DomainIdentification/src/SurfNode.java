package braininferencesurface;
//class representing a surface node (from the .gii file)

import java.util.Vector;


public class SurfNode {
    int nodeID;//0 to N-1
    Vector<Integer> edges; //contains ids of spatially adjacent nodes.
    double[] timeseries;
    boolean hasData = true;
    
    Vector<Integer> localNeighborhood=null;   
    double localNeighScore = -1;
    boolean isPeak=false;
    
    double score;
    boolean updateScore = false;
    
    public SurfNode(int nodeID)
    {
        this.nodeID = nodeID;
        this.edges = new Vector<Integer>();
    }
    
    public void addLocalNeighbors(Vector<Integer> localNeihgs)
    {
        this.localNeighborhood = new Vector<Integer>();
        for(int i = 0; i < localNeihgs.size(); i++)
            localNeighborhood.add(localNeihgs.get(i));
    }
    
    public boolean hasData()
    {
        int cc = 0;
        for(int i = 0; i < timeseries.length; i++)
            if(timeseries[i] == 0)
                cc++;
        if(cc == timeseries.length)
            return false;
        else return true;
    }
    
    public void removeEdge(Integer edge)
    {
        edges.remove(edge);
    }
        
    public void addEdge(Integer edge)
    {
        if(!edges.contains(edge))
            edges.add(edge);
    }
    
    public void addNetworkTimeseriesFromFile(String[] ts, int dimT)
    {
        this.timeseries = new double[dimT];
        for(int i = 0; i < ts.length; i++)
            this.timeseries[i] = Double.parseDouble(ts[i]);
    }
    
    public void addAreaTimeseriesFromFile(String[] ts, int dimT, int maxLag)
    {
        this.timeseries = new double[dimT-2*maxLag];
        int cc = 0;
        for(int i = maxLag; i < ts.length-maxLag; i++)
        {
            this.timeseries[cc] = Double.parseDouble(ts[i]);
            cc++;
        }
    }
    
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 59 * hash + this.nodeID;
        return hash;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(this.getClass() != o.getClass())
            return false;
        SurfNode n = (SurfNode)o;
        return this.nodeID == n.nodeID;
    }
    
    
}
