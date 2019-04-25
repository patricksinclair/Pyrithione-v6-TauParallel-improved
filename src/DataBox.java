public class DataBox {

    private double[] popSizes;
    private double[][] popDistbs;
    private double[] biofilmEdges;
    private double[][] avgGenotypeDistbs;
    private double[][] genoStDevs;

    public DataBox(double[] popSizes, double[][] popDistbs, double[] biofilmEdges, double[][] avgGenotypeDistbs, double[][] genoStDevs){
        this.popSizes = popSizes;
        this.popDistbs = popDistbs;
        this.biofilmEdges = biofilmEdges;
        this.avgGenotypeDistbs = avgGenotypeDistbs;
        this.genoStDevs = genoStDevs;
    }

    public double[] getPopSizes(){
        return popSizes;
    }

    public double[][] getPopDistbs(){
        return popDistbs;
    }

    public double[] getBiofilmEdges(){
        return biofilmEdges;
    }

    public double[][] getAvgGenotypeDistbs(){
        return avgGenotypeDistbs;
    }

    public double[][] getGenoStDevs(){
        return genoStDevs;
    }
}
