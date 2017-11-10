package braininferencesurface;

import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class Area {

    int id,size;
    double score;
    Vector<Integer> nodes;
    boolean expanded = true;
    
    double[] timeseries;
    Vector<AreaEdge> edges = new Vector<AreaEdge>();
    
    public Area(int id)
    {
        this.id = id;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 79 * hash + this.id;
        return hash;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(this.getClass()!=o.getClass())
            return false;
        Area a = (Area)o;
        return a.id == this.id;
      
    }
    
    
}
