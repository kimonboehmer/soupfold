package experiments;

import algorithms.DP;
import algorithms.PartitionFunction;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;

public class Experiments {
    public static double[][][][] basePairProbabilities(StrandPool sp, int m, int sampleSize){
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        dp.computeMFE();
        int maxStrandLength = 0;
        for (int i = 0; i < sp.getNumStrands(); i++) maxStrandLength = Math.max(maxStrandLength, sp.getStrandLength(i));
        double[][][][] table = new double[sp.getNumStrands()][maxStrandLength][sp.getNumStrands()][maxStrandLength];
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure st = dp.backtrack();
            for (int mm = 0; mm < m; mm++){
                int s = st.getStrandFromPosition(mm);
                for (int i = 0; i < sp.getStrandLength(s); i++){
                    int r = st.getPairedStrand(mm, i);
                    if (r >= 0){
                        table[s][i][st.getStrandFromPosition(r)][st.getPairedIndex(mm, i)]++;
                    }
                }
            }
        }
        for (int s = 0; s < sp.getNumStrands(); s++) for (int i = 0; i < sp.getStrandLength(s); i++) for (int r = 0; r < sp.getNumStrands(); r++) for (int j = 0; j < sp.getStrandLength(r); j++){
            table[s][i][r][j] /= sampleSize;
        }
        return table;
    }
    public static double [][] bpTypes(StrandPool sp, int m, int sampleSize){
        double[][] data = new double[sp.getNumStrands()][sp.getNumStrands()];
        double[][][][] table = basePairProbabilities(sp, m, sampleSize);
        for (int s = 0; s < sp.getNumStrands();s++) for (int r = 0; r < sp.getNumStrands(); r++){
            double sum = 0;
            for (int i = 0; i < sp.getStrandLength(s);i++) for (int j = 0; j < sp.getStrandLength(r);j++){
                sum += table[s][i][r][j];
            }
            if (s == r) sum /= 2;
            data[s][r] = sum;
            System.out.printf("From strand %s to strand %s: %f base pairs in average.\n", s, r, sum);
        }
        return data;
    }
    public static double[][][][] basePairProbabilitiesByPosition(StrandPool sp, int m, int sampleSize){
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        dp.computeMFE();
        int maxStrandLength = 0;
        for (int i = 0; i < sp.getNumStrands(); i++) maxStrandLength = Math.max(maxStrandLength, sp.getStrandLength(i));
        double[][][][] table = new double[m][maxStrandLength][m][maxStrandLength];
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure st = dp.backtrack();
            for (int mm = 0; mm < m; mm++) for (int i = 0; i < sp.getStrandLength(st.getStrandFromPosition(mm)); i++){
                int r = st.getPairedStrand(mm, i);
                if (r >= 0){
                    table[mm][i][r][st.getPairedIndex(mm, i)]++;
                }
            }
        }
        for (int s = 0; s < m; s++) for (int i = 0; i < table[s].length; i++) for (int r = 0; r < m; r++) for (int j = 0; j < table[s][i][r].length; j++){
            table[s][i][r][j] /= sampleSize;
        }
        return table;
    }
    public static double interactionProbability(StrandPool sp, int m, int sampleSize){
        double sum = 0;
        int norm = 0;
        double[][][][] table = basePairProbabilitiesByPosition(sp, m, sampleSize);
        for (int s = 0; s < m; s++) for (int i = 0; i < table[s].length; i++){
            double countAll = 0;
            double countInteriorBPs = 0;
            for (int r = 0; r < m; r++) for (int j = 0; j < table[s][i][r].length; j++){
                countAll += table[s][i][r][j];
            }
            if (countAll > 0){
                for (int j = 0; j < table[s][i][s].length; j++) countInteriorBPs += table[s][i][s][j];
                sum += countInteriorBPs / countAll;
                norm++;
            }
        }
        return 1 - (sum / (double) norm);
    }
    public static double[] expNumOccurencesOfStrands(StrandPool sp, int m, int sampleSize){
        double[] data = new double[sp.getNumStrands()];
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        dp.computeMFE();
        for (int l = 0; l < sampleSize; l++) {
            SecondaryStructure st = dp.backtrack();
            for (int i = 0; i < m; i++) data[st.getStrandFromPosition(i)]++;
        }
        for (int i = 0; i < sp.getNumStrands(); i++){
            data[i] /= sampleSize;
            System.out.println(data[i]);
        }
        return data;
    }
}
