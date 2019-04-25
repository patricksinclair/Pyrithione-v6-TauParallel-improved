import org.apache.commons.math3.distribution.LogNormalDistribution;
import java.util.ArrayList;
import java.util.Random;

public class Microhabitat {

    private double mu = Math.log(7.92016113), sigma = 0.10018864;
    private LogNormalDistribution MIC_distribution = new LogNormalDistribution(mu, sigma);

    Random rand = new Random();

    private int K; //karrying kapacity
    private double c; //conc of antimicrobial
    private double b; //migration rate
    private ArrayList<Double> population; //list of MICs of bacteria

    private boolean surface = false, biofilm_region = false, immigration_zone = false;
    private double threshold_density; //fraction of occupation required for biofilm classification

    private int immigration_counter, replication_counter, detachment_counter;


    public Microhabitat(int K, double c, double b){
        this.K = K;
        this.c = c;
        this.b = b;
        this.threshold_density = 0.2;
        this.population = new ArrayList<>(K);
    }


    public int getN(){return population.size();}
    public double getC(){return c;}
    public boolean isBiofilm_region(){return this.biofilm_region;}
    public boolean isSurface(){return this.surface;}
    public boolean isImmigration_zone(){return this.immigration_zone;}
    public double getThreshold_density(){return threshold_density;}
    public ArrayList<Double> getPopulation(){return population;}

    public void setSurface(boolean surface){this.surface = surface;}
    public void setBiofilm_region(boolean biofilm_region){this.biofilm_region = biofilm_region;}
    public void setImmigration_zone(boolean immigration_zone){this.immigration_zone = immigration_zone;}

    public double fractionFull(){
        //double frac = getN()/(double)K;
        //System.out.println("fraction full: N = "+getN()+"\tK = "+K+"\tfrac = "+frac);
        return getN()/(double)K;
    }




    public boolean atBiofilmThreshold(){
        return fractionFull() >= threshold_density;
    }


    public double migrate_rate(){
        //returns 0.5*b for the microhabitat next to the ship hull, to account for the inability to migrate into the hull
        //also for the microhabitat that's the biofilm edge
        return (surface || immigration_zone) ? 0.5*b : b;
    }

    public double beta(int index){
        return population.get(index);
    }

    public double phi_c(int index){
        double cB = c/beta(index);
        return 1. - (6.*cB*cB)/(5. + cB*cB);
    }

    public double replicationOrDeathRate(int index){
        //TODO this isn't quite right I think. it's allowing for growth above N = K
        return (phi_c(index) > 0.) ? phi_c(index)*(1. - getN()/(double)K) : phi_c(index);
    }


    public double getAvgGenotype(){
        if(getN()==0) return 0.;
        else{
            double sum = 0.;
            for(Double geno : population){
                sum += geno;
            }
            return sum/(double)getN();
        }
    }

    public double getStDevOfGenotype(){
        if(getN() <= 1) return 0.;
        else{
            double mean = getAvgGenotype();
            double sumSq = 0.;

            for(Double geno : population){
                sumSq += (geno-mean)*(geno-mean);
            }
            return Math.sqrt(sumSq/(getN()-1.));
        }
    }


    public void addARandomBacterium_x_N(int n_bacteria){
        for(int i = 0; i < n_bacteria; i++){
            population.add(MIC_distribution.sample());
        }
    }

    public void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++){
            population.add(population.get(index));
        }
    }

    public void addABacterium(double MIC){population.add(MIC);}

    public void removeABacterium(int index){population.remove(index);}



}
