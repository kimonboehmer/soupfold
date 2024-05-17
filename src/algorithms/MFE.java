package algorithms;

import algorithms.DPType;

public class MFE implements DPType {
    @Override
    public double noEffect(){
        return 0;
    }
    @Override
    public double initValue(){
        return 0.0;
    }
    @Override
    public double min(double a, double b) {
        return Math.min(a, b);
    }

    @Override
    public double sum(double a, double b) {
        return a + b;
    }
    @Override
    public double E(){
        return -1;
    }
}
