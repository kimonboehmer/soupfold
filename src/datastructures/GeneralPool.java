package datastructures;
import algorithms.DP;

public class GeneralPool implements StrandPool {

    private final Base[][] strandArray;
    private double[][][][][][] M;
    private final int maxLength;

    public GeneralPool(String[] strandArray){
        this.strandArray = new Base[strandArray.length][];
        int maxLength = 0;
        for (int i = 0; i < strandArray.length; i++){
            int l = strandArray[i].length();
            this.strandArray[i] = new Base[l];
            if (maxLength < l) maxLength = l;
            for (int j = 0; j < l; j++) this.strandArray[i][j] = Base.toBase(strandArray[i].charAt(j));
        }
        this.maxLength = maxLength;
    }
    @Override
    public int getNumStrands() {
        return strandArray.length;
    }

    @Override
    public double getM(int m, int s, int i, int r, int j, boolean c) {
        return M[m][s][i][r][j][c? 1: 0];
    }

    @Override
    public void setM(int m, int s, int i, int r, int j, boolean c, double value) {
        M[m][s][i][r][j][c? 1: 0] = value;
    }

    @Override
    public Base getBase(int s, int pos) {
        return strandArray[s][pos];
    }

    @Override
    public int getStrandLength(int s) {
        return strandArray[s].length;
    }

    @Override
    public void initializeTable(int m, int theta, double initValue) {
        int p = getNumStrands();
        M = new double[m+1][p][maxLength][p][maxLength][2];
        for (int s = 0; s < p; s++) {
            for (int i = 0; i < getStrandLength(s); i++) for (int r = 0; r < p; r++) for (int j = 0; j < getStrandLength(r); j++) {
                        for (int mm = 1; mm <= m; mm++) {
                            M[mm][s][i][r][j][0] = DP.NOT_SET;
                            M[mm][s][i][r][j][1] = DP.NOT_SET;
                        }
                        M[0][s][i][r][j][0] = initValue;
                    }
            for (int i = 0; i < maxLength; i++) for (int t = i; t < Math.min(getStrandLength(s), i + theta + 1); t++){
                M[1][s][i][s][t][0] = initValue;
            }
        }
    }

    @Override
    public String toString(int strand) {
        char[] c = new char[getStrandLength(strand)];
        for (int i = 0; i < getStrandLength(strand); i++) c[i] = strandArray[strand][i].toChar();
        return new String(c);
    }
}
