package braininferencesurfacenet;

import java.util.Vector;

/**
 *
 * @author Foudalis
 */
public class AreaEdge {
    
    double corr;
    double weight;
    int lag;
    int from,to;
    boolean undirected=false;
    
    public AreaEdge()
    {}
    
    public AreaEdge(int from, int to)
    {
        this.from = from;
        this.to = to;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(o.getClass()!=this.getClass())
            return false;
        AreaEdge e = (AreaEdge)o;
        return (this.from == e.from && this.to == e.to);
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 13 * hash + this.from+this.to;        
        return hash;
    }
    
}
