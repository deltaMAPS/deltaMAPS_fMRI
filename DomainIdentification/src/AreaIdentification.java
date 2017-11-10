package braininferencesurface;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;

public class AreaIdentification
{
    NodeScoreComparator nodescoreComp = new NodeScoreComparator();
    AreaScoreComp areaScoreComp = new AreaScoreComp();
 
    double[][] timeseries;
    Vector<SurfNode> surfNodes;
    
    Vector<Area> areas;
    Vector<Integer> peaks;
    
    double threshold;
    int dimT;
    int expandFactor = 6;
    
    public AreaIdentification(Vector<SurfNode> surfNodes, double threshold)
    {
        System.out.println("Area Identification Started.");
        this.surfNodes = surfNodes;
        this.threshold = threshold;
        this.dimT = surfNodes.firstElement().timeseries.length;
        int nNodes = surfNodes.size();
        //we do that just to save on speed 
        this.timeseries = new double[nNodes][dimT];
        for(int i = 0; i < nNodes; i++)
        {
            SurfNode node = surfNodes.get(i);
            timeseries[i]=node.timeseries;
        }
        System.out.println("Loading peaks and constructing areas.");
        System.out.println("Expand factor (K): "+this.expandFactor);
        loadPeaks();
        //finished loading/importing data
        //first step. to speed things up, precalculate area scores
        for(int i = 0; i < areas.size(); i++)
            areas.get(i).score = getScore(areas.get(i).nodes);
 
        while(true)
        {
            System.out.println("Merging starts, current areas: "+areas.size());
            boolean merged = mergeAreas();
            System.out.println("Merging stops, current areas: "+areas.size());
            System.out.println("Expanding areas");
            boolean expanded = expandAreas();
            if(!merged && !expanded)
                break;
        }
        printAreaMap();
        exportAreas();
        printAreaSizeDistribution();
        System.out.println("Finished, final areas:"+areas.size());
    }
    
        
    //search for overlaping or adjacent areas and try to merge them
    private boolean mergeAreas()
    {
        boolean merged = false;
        while(true) //merge as many areas as you can.
        {
            //(1) for each area get it's overlapping or adjacent neighbors.
            //(2) keep the best pair of areas to merge (best pair of areas == largest score OR max overlap)
            double bestScore = -1.0;
            int areaMergeAID = -1, areaMergeBID = -1;
            for(int i = 0; i < areas.size(); i++)
            {
                Area areaA = areas.get(i);
                Vector<Integer> geoConnAreas = getGeoConnectedAreas(areaA);
                for(int j = 0; j < geoConnAreas.size(); j++)
                {
                    int areaBID = geoConnAreas.get(j);
                    if(areaA.id > areaBID) //simple trick if areaA overlaps with areaB then the oposite will hold just get the score for one pair
                    {
                        Area areaB = new Area(areaBID);
                        areaB = areas.get(areas.indexOf(areaB));
                        if(areaA.expanded || areaB.expanded)
                        {
                            double score = getMergedScore(areaA, areaB);
                        
                            if(score > bestScore)
                            {//then that is a better pair of areas to merge
                                areaMergeAID = areaA.id;
                                areaMergeBID = areaB.id;
                                bestScore = score;
                            }
                        }
                    }
                }
            }//end now you should have the best area to merge
            if(bestScore > threshold)
            {   //then I can merge the two areas.
                //remopve them from the areas.
                merged = true;
                Area areaToMergeA = new Area(areaMergeAID);
                areaToMergeA = areas.remove(areas.indexOf(areaToMergeA));
                Area areaToMergeB = new Area(areaMergeBID);
                areaToMergeB = areas.remove(areas.indexOf(areaToMergeB));
                //construct the new area
                Area mergedArea = new Area(areaToMergeA.id);//arbitrary id, dont care.
                mergedArea.nodes = new Vector<Integer>();
                mergedArea.nodes.addAll(getUnion(areaToMergeA, areaToMergeB));
                mergedArea.score = bestScore;
                mergedArea.size = mergedArea.nodes.size();
                areas.add(mergedArea);
            }
            else//stop merging
                break;
        }//finished merging (end while)
        return merged;
    }
    
    private double getMergedScore(Area alpha, Area beta)
    {
        //(1) get the union of the grid cells.
        Vector<Integer> union = getUnion(alpha, beta);
        double score = 0.0;
        int nCells = union.size();
        for(int i = 0; i < nCells; i++)
        {
            double[] tsFrom = timeseries[union.get(i)];
            for(int j = (i+1); j < nCells; j++)
            {
                double[] tsTo = timeseries[union.get(j)];
                score+=pearsonCorrel(tsFrom, tsTo);
            }
        }
        double denom = (double)(nCells*(nCells-1)/2);
        return score/denom;
    }

    private Vector<Integer> getUnion(Area alpha, Area beta)
    {
        HashSet<Integer> union = new HashSet<Integer>();
        union.addAll(alpha.nodes);
        union.addAll(beta.nodes);
        Vector<Integer> toReturn = new Vector<Integer>();
        toReturn.addAll(union);
        return toReturn;
    }
    private Vector<Integer> getGeoConnectedAreas(Area anArea)
    {
        Vector<Integer> geoConnectedAreas = new Vector<Integer>();
        
        //step 1. construct a set that will have this area's grid cels as well as their
        //imediate border grid cells.
        HashSet<Integer> extendedAreaCells = new HashSet<Integer>();
        extendedAreaCells.addAll(anArea.nodes);
        //and add the extended neighborhood here
        int nCells = anArea.nodes.size();
        for(int i = 0; i < nCells; i++)
        {
            SurfNode node = surfNodes.get(anArea.nodes.get(i));
            extendedAreaCells.addAll(node.edges);
        }
        //now search for overlaps
        int nAreas = areas.size();
        for(int i = 0; i < nAreas; i++)
        {
            Area alpha = areas.get(i);
            if(alpha.id!=anArea.id)
            {
                HashSet<Integer> union = new HashSet<Integer>();
                union.addAll(extendedAreaCells);
                union.addAll(alpha.nodes);
                if(union.size() < (alpha.nodes.size()+extendedAreaCells.size()))
                    geoConnectedAreas.add(alpha.id);
            }
        }
        return geoConnectedAreas;
    }
    
    private boolean expandAreas()
    {        
        
        boolean startMerging = false;
        boolean expanded = false;
        while(!startMerging)
        {
            expanded = false;
            for(int i = 0; i < areas.size(); i++)
                areas.get(i).score = getScore(areas.get(i).nodes);
            //sort the areas by score in sescending order.
            java.util.Collections.sort(areas,areaScoreComp);
            //try to expand areas, starting from the area with the highest score                
            int nAreas = areas.size(); 
            Area expandedArea = null;
            for(int i = 0; i < nAreas; i++)
            {
                Area area = areas.get(i);
                System.out.println("Expanding area: "+area.id+" score: "+area.score);
                expandedArea = expandArea(area);                
                if(expandedArea.nodes.size() > area.nodes.size() ) //then the area has been expanded
                {
                    System.out.println("Expanded area with score: "+expandedArea.score);
                    System.out.println("Expanded area size: "+expandedArea.nodes.size());
                    //remove the old area and add the new one
                    expanded = true;
                    expandedArea.expanded = true;
                    expandedArea.size = expandedArea.nodes.size();
                    int index = areas.indexOf(area);
                    areas.remove(index);
                    areas.add(index, expandedArea);                    
                    Vector<Integer> neighIDs = getGeoConnectedAreas(expandedArea);
                    System.out.println("Start merging?");
                    for(int j = 0; j < neighIDs.size(); j++)
                    {
                        Area foo = new Area(neighIDs.get(j));
                        foo = areas.get(areas.indexOf(foo));
                        double score = getMergedScore(expandedArea, foo);
                        if(score>threshold)
                        {
                            startMerging=true;
                            System.out.println("Ya.");
                            break;   
                        }
                        else
                            System.out.println("Nope");
                    }
                }          
                if(startMerging)
                    break;
            }
            if(!expanded)
                break;
        }
        return expanded;
    }
    


    
    private Area expandArea(Area area)
    {
        //initialize a set that will hold all grid cells in this area
        HashSet<Integer> inArea = new HashSet<Integer>(); 
        //and add all that are in already
        inArea.addAll(area.nodes);
        //initialize the area
        Area expandedArea = new Area(area.id);
        
        
        expandedArea.nodes = new Vector<Integer>();        
        //add the ones that are inside already
        expandedArea.nodes.addAll(area.nodes);
        expandedArea.score = area.score;//getScore(expandedArea.cells);
        //find the geographically connected neighbors
        Vector<SurfNode> neighboringNodes = new Vector<SurfNode>(); //neighboring cell id's are adj. matrix id's
        for(int i = 0; i < expandedArea.nodes.size(); i++)
        {
            int nodeID = expandedArea.nodes.get(i);
            SurfNode node = surfNodes.get(nodeID);//that works because nodes are sorted
            for(int j = 0; j < node.edges.size(); j++)
            {
                Integer neighID = node.edges.get(j);
                if(!inArea.contains(neighID))
                {
                    SurfNode n = new SurfNode(neighID);
                    if(!neighboringNodes.contains(n))
                        neighboringNodes.add(n);
                }
            }
        }
        
        int nodesEntered = 0;
        //ready to start expanding area
        //Stop criterion I:  best node score is less than the threshold
        //Stop criterion II: when we do not have any neighboring cells left
        //Stop criterion III: when you've added enough nodes.
        while(neighboringNodes.size() > 0) //Stop criterion II check here.
        {
            if(nodesEntered == expandFactor)//Stop criterion III
                break;//added enough nodes.
            //calculate the score of all the nodes
            for(int i = 0; i < neighboringNodes.size(); i++)            
                if(!neighboringNodes.get(i).updateScore)//check if I already have the current node's score
                    neighboringNodes.get(i).score = getScore(neighboringNodes.get(i).nodeID,expandedArea);            
            //sort the nodes by score
            java.util.Collections.sort(neighboringNodes,nodescoreComp);
            //and get the node with the highest score
            SurfNode bestCandidate = neighboringNodes.remove(0);
            
            
            if(bestCandidate.score > threshold)//Stop criterion I check here
            {
                double areaPairs = expandedArea.nodes.size()*(expandedArea.nodes.size()-1)/2;
                double currentAreaScore = areaPairs*expandedArea.score;
                double bcScore = bestCandidate.score*expandedArea.nodes.size();
                expandedArea.score = (bcScore+currentAreaScore)/ (double)(expandedArea.nodes.size()*(expandedArea.nodes.size()+1)/2);
                
                //ok he passed add him to the area
                expandedArea.nodes.add(bestCandidate.nodeID);
                //make him unavailable 
                inArea.add(bestCandidate.nodeID);
                //upadte the number of nodes in area
                nodesEntered++;
                //update the score of nodes
                for(int x = 0; x < neighboringNodes.size(); x++)                                    
                {
                    SurfNode toUpdate = neighboringNodes.get(x);
                    neighboringNodes.get(x).score = updateScore(toUpdate.nodeID,bestCandidate.nodeID, expandedArea.nodes.size()-1,toUpdate.score);            
                    neighboringNodes.get(x).updateScore=true;
                }
                //get his geographically connected neighbors and add them as candidates for expansion
                bestCandidate = surfNodes.get(surfNodes.indexOf(bestCandidate));
                for(int i = 0; i < bestCandidate.edges.size(); i++)
                {
                    Integer neighID = bestCandidate.edges.get(i);
                    if(!inArea.contains(neighID))
                    {
                        SurfNode n = new SurfNode(neighID);
                        if(!neighboringNodes.contains(n))
                            neighboringNodes.add(n);
                    }
                }
            }
            else break;
        }
        System.out.println("Added an area with size: "+expandedArea.nodes.size());
        return expandedArea;//return the area
    }

    
    private double updateScore(int neighbor, int bestCandidate,  int areaSize, double currentScore)
    {
        double[] tsFrom = timeseries[neighbor];
        double[] tsTo = timeseries[bestCandidate];
        double corr = pearsonCorrel(tsFrom, tsTo);
        return (areaSize*currentScore +corr)/(double)(areaSize+1);
    }
    
    private double getScore(int from ,Area areaTo)
    {
        double score = 0;
        double[] tsFrom = timeseries[from];
        for(int i = 0; i < areaTo.nodes.size(); i++)
        {
            double tsTo[] = timeseries[areaTo.nodes.get(i)];
            score += pearsonCorrel(tsFrom, tsTo);
        }
        return score/=(double)areaTo.nodes.size();
    }

    
    private double getScore(Vector<Integer> areaCells)
    {
        int nCells = areaCells.size();
        double score = 0;
        double denom = nCells*(nCells-1)/2;
        
        for(int i = 0; i < nCells; i++)
        {
            double[] ts1 = timeseries[areaCells.get(i)];
            for(int j = (i+1); j < nCells; j++)
            {
                double[] ts2 = timeseries[areaCells.get(j)];
                score+=pearsonCorrel(ts1, ts2);
            }
        }
        return (score/denom);
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
    

    //DEBUG: DONE.
    private void loadPeaks()
    {
        areas = new Vector<Area>();
        peaks = new Vector<Integer>();
        int areaID = 0;
        //search throughthe surfNodes find the peaks and add 
        //an area including the peak and the local neighborhood
        int nNodes = surfNodes.size();
        for(int i = 0; i < nNodes; i++)
        {
            SurfNode node = surfNodes.get(i);
            if(node.isPeak)
            {
                peaks.add(node.nodeID);
                Area newArea = new Area(areaID);
                newArea.nodes = new Vector<Integer>();
                for(int j = 0; j < node.localNeighborhood.size(); j++)
                    newArea.nodes.add(node.localNeighborhood.get(j));
                newArea.size = newArea.nodes.size();
                areas.add(newArea);
                areaID++;
            }
        }
        System.out.println("Initial areas constructed: "+areas.size());
    }

    private void printAreaSizeDistribution()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("areaSizeDistribution.txt"));
            for(int i = 0; i < areas.size(); i++)
                out.println(areas.get(i).nodes.size());
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    private void printAreaMap()
    {
        double[] map = new double[surfNodes.size()];
        for(int i = 0; i < areas.size(); i++)
        {
            double areaColor = Math.random();
            Area alpha = areas.get(i);
            for(int j = 0; j < alpha.nodes.size(); j++)
            {
                int cellid = alpha.nodes.get(j);
                if(map[cellid] == 0)
                    map[cellid] = areaColor;
                else
                    map[cellid] = -1.0;
            }
        }
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("areasOut.txt"));
            for(int i = 0; i < map.length; i++)
                out.println(map[i]);
            out.close();
        }
        catch(IOException ioe)
        {}
    }
    
    private void exportAreas()
    {
        try
        {
            PrintWriter out = new PrintWriter(new FileWriter("finalAreas.txt"));
            for(int i = 0; i < areas.size(); i++)
            {
                Area alpha = areas.get(i);
                out.println(alpha.id+"\t"+alpha.size+"\t"+alpha.score);
                for(int j = 0; j < alpha.nodes.size(); j++)
                    out.print(alpha.nodes.get(j)+"\t");
                out.println();
            }
            out.close();
        }
        catch(IOException ioe)
        {}
    }

}

class NodeScoreComparator implements Comparator<SurfNode>
{
    @Override
    public int compare(SurfNode n1, SurfNode n2)
    {
        if(n1.score < n2.score)
            return 1;
        else if(n1.score > n2.score)
            return -1;
        else return 0;
                                                
    }
}

class AreaScoreComp implements Comparator<Area>
{
    @Override
    public int compare(Area n1, Area n2)
    {
        if(n1.score < n2.score)
            return 1;
        else if(n1.score > n2.score)
            return -1;
        else return 0;
                                                
    }
}
