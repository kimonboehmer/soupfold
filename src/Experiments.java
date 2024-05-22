
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
                    int r = st.getPairedStrand(s, i);
                    if (r >= 0){
                        table[s][i][r][st.getPairedIndex(s, i)]++;
                    }
                }
            }
        }
        for (int s = 0; s < sp.getNumStrands(); s++) for (int i = 0; i < sp.getStrandLength(s); i++) for (int r = 0; r < sp.getNumStrands(); r++) for (int j = 0; j < sp.getStrandLength(r); j++){
            table[s][i][r][j] /= sampleSize;
        }
        return table;
    }
}
