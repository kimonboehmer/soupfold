public interface DPType {
    double min(double a, double b);
    double sum(double a, double b);
    double noEffect();
    double initValue();
    double E();
    double INFTY = 100000;
    boolean btChoose(double a);
    void btInit(double a);
}
