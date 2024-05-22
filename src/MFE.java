public class MFE implements DPType{
    double targetVal;
    @Override
    public double forbidden(){
        return INFTY;
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

    @Override
    public boolean btChoose(double a) {
        return a == targetVal;
    }

    @Override
    public void btInit(double a) {
        targetVal = a;
    }

    @Override
    public double strandPenalty(int length) {
        return 0;
    }


}
