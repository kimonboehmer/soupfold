
import experiments.Helper;

import java.io.IOException;

import static experiments.FinalExperiments.doExperiment;

public class Main {
    /**
     * call doExperiment() with the number of the experiment to reproduce.
     * 0: all experiments below
     * 1: MFE structure for homogeneous strand soup
     * 2: homogeneous strand soup MFE comparing all TR patterns
     * 3: MFE depending on sequence length
     * 4: MFE depending on number of interacting strands
     * 5: MFE structure for heterogeneous strand soup
     * 6: exterior homo/hetero base pair probability for one pair of TR patterns, with increasing m
     * 7: main experiment: exterior homo/hetero base pair probability for all pairs of TR patterns
     * 8: speed test: Measures the time the algorithm needs for increasing m.
     * 9: gives a matrix of the pairwise bp probabilities of the "winning" MIS strand soup for m=10
     */
    public static void main(String[] args) throws IOException {
        doExperiment(0);
    }
}
