import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

public class BioSystem {

    Random rand = new Random();

    private int L, K; //length, karrying kapacity
    private double alpha, c_max; //steepness of antimicrobial gradient, max concn
    private Microhabitat[] microhabitats;
    private double timeElapsed;
    private double tau; //timestep used in tau-leaping
    private double immigration_rate = 7.2;
    private double migration_rate = 0.2;
    private double attachment_rate = 2000.;
    private double detachment_rate = 1.38*0.1;
    private double delta_x = 5.;
    private int immigration_index, biofilm_edge_index;


    public BioSystem(int L, int K, double alpha, double c_max, double tau){
        this.L = L;
        this.K = K;
        this.alpha = alpha;
        this.c_max = c_max;
        this.tau = tau;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;
        this.immigration_index = 0;

        for(int i = 0; i < L; i++) {
            double c_i = c_max*Math.exp(-alpha*i*delta_x);
            microhabitats[i] = new Microhabitat(K, c_i, migration_rate);
        }

        microhabitats[0].setSurface(true);
        microhabitats[0].addARandomBacterium_x_N(25);
    }


    public double getTimeElapsed(){
        return timeElapsed;
    }

    public int getImmigration_index(){return immigration_index;}

    public int getPop_i(int index){
        return microhabitats[index].getN();
    }

    public int getN_i(int index){
        return microhabitats[index].getN();
    }

    public int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    public double[] getConcentrationProfile(){
        double[] c_vals = new double[L];
        for(int i = 0; i < L; i++){
            c_vals[i] = microhabitats[i].getC();
        }
        return c_vals;
    }

    public boolean[] getBiofilmProfile(){
        boolean[] bf_vals = new boolean[L];
        for(int i = 0; i < L; i++){
            bf_vals[i] = microhabitats[i].atBiofilmThreshold();
        }
        return bf_vals;
    }


    public int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = 0; i < L; i++) {
            //changed this so it doesn't have a break statement
            if(microhabitats[i].isBiofilm_region()) {
                edgeIndex = i;
                //System.out.println(edgeIndex);
                //break;
            }
        }
        return edgeIndex;
    }


    public int getBiofilmSize(){
        return getBiofilmEdge()+1;
    }


    public double[] getPopulationDistribution(){
        double[] popSizes = new double[L];
        for(int i = 0; i < L; i++) {
            popSizes[i] = microhabitats[i].getN();
        }
        return popSizes;
    }


    public double[] getAvgGenotypeDistribution(){
        double[] avgGenos = new double[L];
        for(int i = 0; i < L; i++) {
            avgGenos[i] = microhabitats[i].getAvgGenotype();
        }
        return avgGenos;
    }


    public double[] getStDevOfGenotypeDistribution(){
        double[] genoStDevs = new double[L];
        for(int i = 0; i < L; i++) {
            genoStDevs[i] = microhabitats[i].getStDevOfGenotype();
        }
        return genoStDevs;
    }


    public void immigrate(int mh_index, int n_immigrants){
        microhabitats[mh_index].addARandomBacterium_x_N(n_immigrants);
    }



    public void migrate(Microhabitat[] updated_microhabs, int mh_index, int bac_index){
        //TODO this needs seriously reworked to handle cases of only one microhabitat - done with if statement in update method
        //TODO allow migration from immigration zone into biofilm before immigration zone -> biofilm? -done

        int biof_edge = immigration_index;
        double migrating_bac = updated_microhabs[mh_index].getPopulation().get(bac_index);
        updated_microhabs[mh_index].removeABacterium(bac_index);

        if(microhabitats[mh_index].isSurface()){
            updated_microhabs[mh_index+1].addABacterium(migrating_bac);

        }else if(microhabitats[mh_index].isImmigration_zone()){
            updated_microhabs[mh_index-1].addABacterium(migrating_bac);

        }else{
            if(rand.nextBoolean()){
                updated_microhabs[mh_index+1].addABacterium(migrating_bac);
            }else{
                updated_microhabs[mh_index-1].addABacterium(migrating_bac);
            }
        }
    }


    public void updateBiofilmSize(){
        //once the edge microhabitat is sufficiently populated, this moves the immigration edge
        //microhabitat to the next microhabitat along
        int i_counter = 0; //no idea why i need to use i_counter instead of immigration_index=getBiofilmEdge()+1
        int i_lim = Math.min(L, immigration_index+5);
        for(int i = 0; i < L; i++){

            microhabitats[i].setImmigration_zone(false);

            if(microhabitats[i].atBiofilmThreshold()){

                microhabitats[i].setBiofilm_region(true);
                //biofilm_edge_index = i;
                i_counter = i+1;

            }
        }
        //maybe incorporate a limit on this in case it ever reaches L
        immigration_index = Math.max(immigration_index, i_counter); //this stops the immigration index regressing if a microhabitat depletes
        microhabitats[immigration_index].setImmigration_zone(true);
    }



    public void performAction(){

        Microhabitat[] updated_microhabs = microhabitats.clone(); //the new system that will replace the current one at each iteration

        double tau_step = tau;
        int biofilm_thickness = getBiofilmSize(); //the +1 includes the immigration microhab.
        //int biofilm_edge = getBiofilmEdge(); //this is the furthest index of biofilm. not including immig index
        int[][] replication_allocations;
        int[][] death_allocations;
        int[][] migration_allocations;
        int[] detachment_allocations;
        int[] original_popsizes;
        int n_immigrants;


        //whileloop to calculate all the reactions to take place. used to carry out modified tau-leaping
        whileloop:
        while(true){

            replication_allocations = new int[biofilm_thickness][];
            death_allocations = new int[biofilm_thickness][];
            migration_allocations = new int[biofilm_thickness][];
            original_popsizes = new int[biofilm_thickness];
            detachment_allocations = new int[microhabitats[immigration_index].getN()];


            for(int mh_index = 0; mh_index < biofilm_thickness; mh_index++){
                int mh_pop = microhabitats[mh_index].getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];
                int[] n_migrations = new int[mh_pop];


                for(int bac_index = 0; bac_index < mh_pop; bac_index++){


                    ////////// MIGRATIONS //////////////////////
                    n_migrations[bac_index] = new PoissonDistribution(microhabitats[mh_index].migrate_rate()*tau_step).sample();

                    if(n_migrations[bac_index] > 1){
                        tau_step /= 2.;
                        continue whileloop;
                    }
                    ////////////////////////////////////////////



                    ///////////// DETACHMENTS /////////////////////////
                    if(mh_index == immigration_index){
                        detachment_allocations[bac_index] = new PoissonDistribution(detachment_rate*tau_step).sample();
                        if(detachment_allocations[bac_index] > 1){
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a bacteria is detaching then it can't migrate
                        if(detachment_allocations[bac_index] != 0){
                            n_migrations[bac_index] = 0;
                        }
                    }
                    ////////////////////////////////////////////////////////



                    ////////////////// REPLICATIONS AND DEATHS ///////////////////////////
                    double g_or_d_rate = microhabitats[mh_index].replicationOrDeathRate(bac_index);

                    if(g_or_d_rate == 0.){

                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = 0;

                    }else if(g_or_d_rate > 0){

                        n_replications[bac_index] = new PoissonDistribution(g_or_d_rate*tau_step).sample();
                        n_deaths[bac_index] = 0;

                    }else{
                        n_replications[bac_index] = 0;
                        n_deaths[bac_index] = new PoissonDistribution(Math.abs(g_or_d_rate)*tau_step).sample();

                        if(n_deaths[bac_index] > 1){
                            tau_step /= 2.;
                            continue whileloop;
                        }
                        //if a death is occurring, then that bacteria can't migrate or detach
                        if(n_deaths[bac_index] !=0) {
                            n_migrations[bac_index] = 0;
                            if(mh_index == immigration_index) detachment_allocations[bac_index] = 0;
                        }
                    }
                    /////////////////////////////////////////////////////////////////////////
                }

                replication_allocations[mh_index] = n_replications;
                death_allocations[mh_index] = n_deaths;
                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = microhabitats[mh_index].getN();
            }

            n_immigrants = new PoissonDistribution(immigration_rate*tau_step).sample();
            break whileloop;
        }


        for(int mh_index = 0; mh_index < biofilm_thickness; mh_index++){

            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(death_allocations[mh_index][bac_index] != 0) updated_microhabs[mh_index].removeABacterium(bac_index);

                else{
                    updated_microhabs[mh_index].replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);

                    if(getBiofilmSize() > 1){
                        //this statement ensures migration can only occur between biofilm microhabitats
                        if(migration_allocations[mh_index][bac_index] != 0) migrate(updated_microhabs, mh_index, bac_index);
                    }

                    if(mh_index == immigration_index){
                        if(detachment_allocations[bac_index] != 0) updated_microhabs[mh_index].removeABacterium(bac_index);
                    }
                }
            }
        }

        microhabitats = updated_microhabs;
        immigrate(immigration_index, n_immigrants);
        updateBiofilmSize();
        timeElapsed += tau_step;
    }







    public static void tester(){

        double duration = 2000.;
        int L = 500;
        int K = 500;
        double c_max = 10.;
        double alpha = 1e-4;
        double tau = 0.01;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);
        System.out.println(Arrays.toString(bs.getConcentrationProfile())+"\n");

        while(bs.getTimeElapsed() <= duration) {

            if(bs.timeElapsed%10. < 0.01) {

                System.out.println("-----------------------------------------------------------------------------------------------");

                String output = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop: %d \tbf_edge fracfull: %.3f \timmig index: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bs.getBiofilmEdge()), bs.microhabitats[bs.getBiofilmEdge()].fractionFull(), bs.getImmigration_index());

                String output2 = String.format("time elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d",
                        bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge());

                System.out.println(output);
                System.out.println("\nN distb");
                System.out.println(Arrays.toString(bs.getPopulationDistribution()) + "\n");
                System.out.println("Biofilm y/n?");
                System.out.println(Arrays.toString(bs.getBiofilmProfile()) + "\n");
                System.out.println("-----------------------------------------------------------------------------------------------");
            }
            bs.performAction();

        }

    }


    public static void getPopDistbInfo(){
        //method to get info on population distbs
        //get popsize over time
        //pop distb over time
        //biofilm edge over time
        //avg genotype distb over time

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 10, nMeasurements = 20;
        double duration = 100., interval = duration/nMeasurements;

        String popSizeFilename = "pyrithione-testing-pop_size-t=" + String.valueOf(duration);
        String popDistbFilename = "pyrithione-testing-pop_distb-t=" + String.valueOf(duration);
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t=" + String.valueOf(duration);
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t=" + String.valueOf(duration);
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t=" + String.valueOf(duration);
        String counterDistbsFilename = "pyrithione-testing-counterDistb-t=" + String.valueOf(duration);

        String[] counterHeaders = {"immigration", "migrationIn", "migrationOut", "replication", "death"};
        int nCounters = counterHeaders.length;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];
        double[][] allCounters = new double[nReps][];


        for(int r = 0; r < nReps; r++) {

            BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);
            //double[] cVals = bs.getCVals();
            //System.out.println(Arrays.toString(cVals));
            //System.out.println(cVals[476]);

            boolean alreadyRecorded = false;
            int timerCounter = 0;

            double[] popSizes = new double[nMeasurements + 1];
            double[][] popDistbs = new double[nMeasurements + 1][];
            double[] biofilmEdges = new double[nMeasurements + 1];
            double[][] avgGenotypeDistbs = new double[nMeasurements + 1][];
            double[][] genoStDevs = new double[nMeasurements + 1][];


            while(bs.timeElapsed <= duration + 0.02*interval) {

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded) {

                    int bf_edge_i = bs.getBiofilmEdge();

                    String output = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d \tbf_edge pop %d \tbf_edge fracfull: %.3f",
                            r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bf_edge_i), bs.microhabitats[bf_edge_i].fractionFull());

                    String output2 = String.format("rep: %d \ttime elapsed: %.3f \ttotal N: %d \tbiofilm edge: %d, bf_edge pop: %d",
                            r, bs.getTimeElapsed(), bs.getTotalN(), bs.getBiofilmEdge(), bs.getN_i(bf_edge_i));

                    System.out.println(output2);

                    popSizes[timerCounter] = bs.getTotalN();
                    popDistbs[timerCounter] = bs.getPopulationDistribution();
                    biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                    avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                    genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();

                    alreadyRecorded = true;
                    timerCounter++;
                }
                if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

                bs.performAction();
            }

            allPopSizes[r] = popSizes;
            allPopDistbs[r] = popDistbs;
            allBiofilmEdges[r] = biofilmEdges;
            allAvgGenotypeDistbs[r] = avgGenotypeDistbs;
            allGenoStDevs[r] = genoStDevs;
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        System.out.println("results written to file");
    }



    public static double[][] getAvgGenoDistbs(int i){

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nMeasurements = 20;
        double duration = 100., interval = duration/nMeasurements;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        boolean alreadyRecorded = false;
        int timerCounter = 0;

        double[][] avgGenotypeDistbs = new double[nMeasurements + 1][];

        while(bs.timeElapsed <= duration + 0.02*interval) {

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.01*interval) && !alreadyRecorded) {

                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed());
                avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                alreadyRecorded = true;
                timerCounter++;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        return avgGenotypeDistbs;
    }





    public static DataBox getAllData(int i){

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nMeasurements = 40;
        double duration = 200., interval = duration/nMeasurements;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);

        boolean alreadyRecorded = false;
        int timerCounter = 0;

        double[] popSizes = new double[nMeasurements+1];
        double[][] popDistbs = new double[nMeasurements+1][];
        double[] biofilmEdges = new double[nMeasurements+1];
        double[][] avgGenotypeDistbs = new double[nMeasurements+1][];
        double[][] genoStDevs = new double[nMeasurements+1][];

        while(bs.timeElapsed <= duration+0.02*interval){

            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmSize()*K;
                int total_N = bs.getTotalN();

                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());

                popSizes[timerCounter] = total_N;
                popDistbs[timerCounter] = bs.getPopulationDistribution();
                biofilmEdges[timerCounter] = bs.getBiofilmEdge();
                avgGenotypeDistbs[timerCounter] = bs.getAvgGenotypeDistribution();
                genoStDevs[timerCounter] = bs.getStDevOfGenotypeDistribution();

                alreadyRecorded = true;
                timerCounter++;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        return new DataBox(popSizes, popDistbs, biofilmEdges, avgGenotypeDistbs, genoStDevs);
    }



    public static int getThicknessReachedAfterATime(double duration, int i){
        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        BioSystem bs = new BioSystem(L, K, alpha, c_max, tau);
        int nUpdates = 20;
        double interval = duration/nUpdates;
        boolean alreadyRecorded = false;

        while(bs.timeElapsed <= duration){


            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getBiofilmSize()*K;
                int total_N = bs.getTotalN();
                System.out.println("rep : "+i+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge());

                alreadyRecorded = true;
            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;



            bs.performAction();
        }

        return bs.getBiofilmEdge();
    }


    public static void getBiofilmThicknessHistoInParallel(int nReps){
        long startTime = System.currentTimeMillis();

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        double duration = 1680.; //10 week duration


        int[] mh_index_reached = new int[nReps];
        String index_reached_filename = "pyrithione-testing-mh_index_reached_histo-t="+String.valueOf(duration)+"parallel";

        IntStream.range(0, nReps).parallel().forEach(i -> mh_index_reached[i] = BioSystem.getThicknessReachedAfterATime(duration, i));

        Toolbox.writeHistoArrayToFile(index_reached_filename, mh_index_reached);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }


    public static void getInfoInParallel(){

        long startTime = System.currentTimeMillis();

        int K = 500, L = 500;
        double c_max = 10., alpha = 0.01, tau = 0.01;

        int nReps = 16;
        double duration = 200.;

        double[][] allPopSizes = new double[nReps][];
        double[][][] allPopDistbs = new double[nReps][][];
        double[][] allBiofilmEdges = new double[nReps][];
        double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        double[][][] allGenoStDevs = new double[nReps][][];

        String popSizeFilename = "pyrithione-testing-pop_size-t="+String.valueOf(duration)+"-parallel";
        String popDistbFilename = "pyrithione-testing-pop_distb-t="+String.valueOf(duration)+"-parallel";
        String biofilmEdgeFilename = "pyrithione-testing-biofilm_edge-t="+String.valueOf(duration)+"-parallel";
        String avgGenotypeDistbFilename = "pyrithione-testing-avgGenoDistb-t="+String.valueOf(duration)+"-parallel";
        String genoStDevDistbFilename = "pyrithione-testing-genoStDevDistb-t="+String.valueOf(duration)+"-parallel";

        //double[][][] allAvgGenotypeDistbs = new double[nReps][][];
        DataBox[] dataBoxes = new DataBox[nReps];

        IntStream.range(0, nReps).parallel().forEach(i -> dataBoxes[i] = BioSystem.getAllData(i));

        for(int j = 0; j < dataBoxes.length; j++){
            allPopSizes[j] = dataBoxes[j].getPopSizes();
            allPopDistbs[j] = dataBoxes[j].getPopDistbs();
            allBiofilmEdges[j] = dataBoxes[j].getBiofilmEdges();
            allAvgGenotypeDistbs[j] = dataBoxes[j].getAvgGenotypeDistbs();
            allGenoStDevs[j] = dataBoxes[j].getGenoStDevs();
        }

        double[] processedPopSizes = Toolbox.averagedResults(allPopSizes);
        double[][] processedPopDistbs = Toolbox.averagedResults(allPopDistbs);
        double[] processedBiofilmEdges = Toolbox.averagedResults(allBiofilmEdges);
        double[][] processedAvgGenotypeDistbs = Toolbox.averagedResults(allAvgGenotypeDistbs);
        double[][] processedGenoStDevs = Toolbox.averagedResults(allGenoStDevs);

        Toolbox.writeAveragedArrayToFile(popSizeFilename, processedPopSizes);
        Toolbox.writeAveragedDistbsToFile(popDistbFilename, processedPopDistbs);
        Toolbox.writeAveragedArrayToFile(biofilmEdgeFilename, processedBiofilmEdges);
        Toolbox.writeAveragedDistbsToFile(avgGenotypeDistbFilename, processedAvgGenotypeDistbs);
        Toolbox.writeAveragedDistbsToFile(genoStDevDistbFilename, processedGenoStDevs);

        long finishTime = System.currentTimeMillis();

        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);

        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }







    public static double getCValWithOffset(int index, double maxC, double alpha, int L){
        //this calculates i* for the gradient profile offset, moves so the final concn is maxC, and it decreases with 1/e
        //or something like that
        //then calculates the corresponding concentration in that microhabitat

        double offset =  (L-1.) - Math.log(maxC+1.)/alpha;
        return (index >= offset) ? Math.exp(alpha*(index - offset)) - 1. : 0.;
    }
}
