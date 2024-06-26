package algorithms;

import algorithms.DPType;

import java.util.Random;

public class PartitionFunction implements DPType {
    private double randomVal;
    private final Random r;
    double bpContrib;
    final double BOLTZMANN_CONSTANT = 0.001987;

    /**
     * @param temperature in Kelvin
     */
    public PartitionFunction(double temperature){
        this.bpContrib = Math.exp(1.0 / (temperature * BOLTZMANN_CONSTANT));
        r = new Random();
    }
    public double forbidden(){
        return 0;
    }
    public double initValue() {
        return 1.0;
    }
    @Override
    public double min(double a, double b) {
        return a + b;
    }

    @Override
    public double sum(double a, double b) {
        return a * b;
    }
    @Override
    public double E(){
        return bpContrib;
    }

    @Override
    public boolean btChoose(double a) {
        randomVal -= a;
        return randomVal <= 0;
    }

    @Override
    public void btInit(double a) {
        randomVal = r.nextDouble(a);
    }

    @Override
    public double strandPenalty(int strand, int length) {
        return Math.pow(1, length) ;
    }

}
