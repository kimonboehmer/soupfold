package algorithms;

import algorithms.DPType;

public class PartitionFunction implements DPType {
    double bpContrib;
    final double BOLTZMANN_CONSTANT = 0.001987;
    public PartitionFunction(double temperature){
        this.bpContrib = Math.exp(- 1.0 / (temperature * BOLTZMANN_CONSTANT));
    }
    public double noEffect(){
        return 1.0;
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
}
