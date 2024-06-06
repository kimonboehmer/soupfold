package datastructures;

import algorithms.DP;

import java.util.*;

public class TripletPool implements StrandPool {
    private static final byte REPEAT_LENGTH = 3;
    protected static class Triplet {
        int pattern;
        int repeats;
        protected Triplet(int pattern, int repeats){
            this.pattern = pattern;
            this.repeats = repeats;
        }
    }
    int maxRepeats;
    int minRepeats;
    int numPatterns;
    private final Base[][] patternArray;
    private final Triplet[] pool;
    private double[][][][][][] M;
    public TripletPool(Base[] pattern, int mid, int rad, int num){
        pool = new Triplet[num];
        Random rand = new Random(1);
        minRepeats = 100000;
        maxRepeats = 0;
        for (int i = 0; i < num; i++){
            int r = rand.nextInt(mid - rad, mid + rad);
            pool[i] = new Triplet(0, r);
            if (minRepeats > r) minRepeats = r;
            if (maxRepeats < r) maxRepeats = r;
        }
        patternArray = new Base[][]{pattern};
        numPatterns = 1;
    }
    public TripletPool(Base[] pattern, int mid, int rad){
        assert mid > rad;
        pool = new Triplet[2 * rad + 1];
        pool[0] = new Triplet(0, mid);
        for (int i = 1; i <= rad; i++){
            pool[2 * i - 1] = new Triplet(0, mid - i);
            pool[2 * i] = new Triplet(0, mid + i);
        }
        patternArray = new Base[][]{pattern};
        numPatterns = 1;
        maxRepeats = mid + rad;
        minRepeats = mid - rad;
    }
    private int tripletIndex(LinkedList<Base[]> ll,Base[] t){
        int i = 0;
        for (Base[] r : ll) {
            if (r[0] == t[0] && r[1] == t[1] && r[2] == t[2]) return i;
            i++;
        }
        return -1;
    }
    public TripletPool(List<String> strands){
        pool = new Triplet[strands.size()];
        maxRepeats = 0;
        minRepeats = 100000;
        LinkedList<Base[]> patternOccurences = new LinkedList<>();
        int i = 0;
        int val = 0;
        for (String strand : strands){
            char[] strandCharArray = strand.toCharArray();
            Base[] pattern = new Base[REPEAT_LENGTH];
            for (int j = 0; j < REPEAT_LENGTH; j++){
                pattern[j] = Base.toBase(strandCharArray[j]);
            }
            if (tripletIndex(patternOccurences, pattern) == -1){
                val++;
                patternOccurences.add(pattern);
            }
            int reps = Integer.parseInt(strand.substring(REPEAT_LENGTH));
            if (maxRepeats < reps) maxRepeats = reps;
            if (minRepeats > reps) minRepeats = reps;
            Triplet t = new Triplet(tripletIndex(patternOccurences, pattern), reps);
            pool[i++] = t;
        }
        numPatterns = val;
        patternArray = new Base[numPatterns][REPEAT_LENGTH];
        i = 0;
        for (Base[] pattern : patternOccurences) patternArray[i++] = pattern;
    }
    @Override
    public void initializeTable(int m, int theta, double initValue){
        M = new double[m+1][numPatterns][maxRepeats * REPEAT_LENGTH][numPatterns][maxRepeats * 3][2];
        for(int si=0;si<numPatterns;si++)for(int ii=0;ii<maxRepeats * REPEAT_LENGTH;ii++)for(int ri=0;ri<numPatterns;ri++)for(int ji=0;ji<maxRepeats * REPEAT_LENGTH;ji++){
            for (int mi = 1; mi < m+1;mi++){ // initialize everything with "uninitialized"
                M[mi][si][ii][ri][ji][0]= DP.NOT_SET;
                M[mi][si][ii][ri][ji][1]= DP.NOT_SET;
            }
            M[0][si][ii][ri][ji][0] = initValue; // if no strands are remaining, there is no contribution
        }
        for (int i = 0; i < REPEAT_LENGTH; i++){
            for (int t = 0; t <= Math.min(theta, maxRepeats * REPEAT_LENGTH - 1); t++){
                for (int p = 0; p < numPatterns; p++){
                    M[1][p][i][p][t][0] = initValue; // single-stranded intervals which do not
                    // respect the minimum base pair span have no contribution
                }
            }
        }
    }
    @Override
    public Base getBase(int s, int pos) {
        if (pos > pool[s].repeats * REPEAT_LENGTH) throw new InputMismatchException("Position greater than strand length!");
        return patternArray[pool[s].pattern][pos % REPEAT_LENGTH];
    }
    @Override
    public double getM(int m, int s, int i, int r, int j, boolean c){
        if (m == 1) return M[1][pool[s].pattern][i % REPEAT_LENGTH][pool[s].pattern][j - i][0];
        return M[m][pool[s].pattern][getStrandLength(s) - i - 1][pool[r].pattern][j][c ? 1 : 0];
    }
    @Override
    public void setM(int m, int s, int i, int r, int j, boolean c, double val){
        if (m == 1) M[1][pool[s].pattern][i % REPEAT_LENGTH][pool[s].pattern][j - i][0] = val;
        M[m][pool[s].pattern][getStrandLength(s) - i - 1][pool[r].pattern][j][c ? 1 : 0] = val;
    }
    @Override
    public int getStrandLength(int s){
        return pool[s].repeats * REPEAT_LENGTH;
    }
    @Override
    public int getNumStrands(){
        return pool.length;
    }
    @Override
    public String toString(int strand){
        char[] triplet = new char[REPEAT_LENGTH * pool[strand].repeats];
        for (int i = 0; i < triplet.length; i++) triplet[i] = patternArray[getPattern(strand)][i % REPEAT_LENGTH].toChar();
        return new String(triplet);
    }
    public int getPattern(int strand){
        return pool[strand].pattern;
    }
}
