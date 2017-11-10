package braininferencesurfacenet;

import java.util.Vector;

public class Area {

    int id;
    int cortex;//cortex = 0 for left cortex, 1 for right cortex
    Vector<Integer> surfVoxels;
    Vector<AreaEdge> edges = new Vector<AreaEdge>();
    
    int size;
    double strength;
    double[] timeseries;
    
    public Area(int id)
    {
        this.id = id;
        surfVoxels = new Vector<Integer>();        
    }
    
    public void setCortex(int cortex)
    {this.cortex = cortex;}
    
    public void addTimeSeries(double[] timeseries)
    {
        this.timeseries = new double[timeseries.length];
        for(int i = 0; i < timeseries.length; i++)
            this.timeseries[i]=timeseries[i];
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(this.getClass()!= o.getClass())
            return false;
        Area a = (Area)o;
        return this.id == a.id;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + this.id;
        return hash;
    }
}
