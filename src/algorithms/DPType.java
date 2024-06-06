package algorithms;

public interface DPType {
    /**
     * @return the minimum of a and b for the MFE, and the sum of a and b for the PF.
     */
    double min(double a, double b);

    /**
     * @return the sum of a and b for the MFE, and the product of a and b for the PF.
     */
    double sum(double a, double b);

    /**
     * @return a value for a forbidden secondary structure.
     */
    double forbidden();

    /**
     * @return a value for an empty secondary structure.
     */
    double initValue();

    /**
     * @return value for the energy contribution of one base pair.
     */
    double E();
    double INFTY = 100000;

    /**
     * @param a value of the current backtracking candidate
     * @return true if we choose this candidate in the backtracking
     */
    boolean btChoose(double a);

    /**
     * @param a total MFE/PF value before a backtracking step.
     */
    void btInit(double a);

    /**
     * @param strand strand number
     * @param length strand length
     * @return a penalty depending on the added strand.
     */
    double strandPenalty(int strand, int length);
}
