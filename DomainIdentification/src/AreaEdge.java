package braininferencesurface;

/**
 *
 * @author ilias.fountalis
 */
public class AreaEdge {
    
    int from, to; 
    double weight;
    double sij;
    double lag;    
    double minLag, maxLag;
    
    public AreaEdge()
    {}
    
    public AreaEdge(int from, int to, double weight, double lag)
    {
        this.from = from;
        this.to = to;
        this.weight = weight;
        this.lag = lag;
    }
    
    

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 11 * hash + this.from;
        hash = 11 * hash + this.to;
        return hash;
    }
    
    @Override
    public boolean equals(Object o)
    {
        if(this.getClass() != o.getClass())
            return false;
        AreaEdge e = (AreaEdge)o;
        return this.from == e.from && this.to == e.to;
    }
}
