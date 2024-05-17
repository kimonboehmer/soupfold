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
    private int minDiff;
    private int maxDiff;
    private final Base[][] patternArray;
    private final Triplet[] pool;
    double[][][][][][][] M;
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
    public TripletPool(List<String> strands){
        pool = new Triplet[strands.size()];
        maxRepeats = 0;
        minRepeats = 100000;
        HashMap<Base[], Integer> hm = new HashMap<>();
        int counter = 0;
        int i = 0;
        for (String strand : strands){
            char[] strandChars = strand.toCharArray();
            Base[] pattern = new Base[REPEAT_LENGTH];
            for (int j = 0; j < REPEAT_LENGTH; j++){
                pattern[j] = Base.toBase(strandChars[j]);
            }
            Integer val = hm.get(pattern);
            if (val == null){
                hm.put(pattern, counter);
                val = counter++;
            }
            int reps = Integer.parseInt(strand.substring(REPEAT_LENGTH));
            if (maxRepeats < reps) maxRepeats = reps;
            if (minRepeats > reps) minRepeats = reps;
            Triplet t = new Triplet(val, reps);
            pool[i++] = t;
        }
        numPatterns = counter;
        patternArray = new Base[numPatterns][REPEAT_LENGTH];
        for (Map.Entry<Base[], Integer> entry : hm.entrySet()){
            patternArray[entry.getValue()] = entry.getKey();
        }
    }
    public void initializeTable(int m, int theta, double initValue, int avg){
        minDiff =  m * (avg - minRepeats*3);
        maxDiff =  m * (maxRepeats*3 - avg);
        int diffRange = 2 * minDiff + 1 + 2* maxDiff;
        M = new double[m+1][numPatterns][maxRepeats * REPEAT_LENGTH][numPatterns][maxRepeats * 3][2][diffRange];
        for(int si=0;si<numPatterns;si++)for(int ii=0;ii<maxRepeats * REPEAT_LENGTH;ii++)for(int ri=0;ri<numPatterns;ri++)for(int ji=0;ji<maxRepeats * REPEAT_LENGTH;ji++){
            for (int mi = 1; mi < m+1;mi++) for (int di = 0; di < diffRange; di++){
                M[mi][si][ii][ri][ji][0][di]= DP.NOT_SET;
                M[mi][si][ii][ri][ji][1][di]= DP.NOT_SET;
            }
            M[0][si][ii][ri][ji][0][2*minDiff] = initValue;
        }
        for (int i = 0; i < maxRepeats * REPEAT_LENGTH; i++){
            for (int t = 0; t < Math.min(theta, maxRepeats * REPEAT_LENGTH - i); t++){
                for (int p = 0; p < numPatterns; p++){
                    M[1][p][i][p][i + t][0][2*minDiff] = initValue;
                }
            }
        }
    }
    public Base getBase(int s, int pos) {
        if (pos > pool[s].repeats * REPEAT_LENGTH) throw new InputMismatchException("Position greater than strand length!");
        return patternArray[pool[s].pattern][pos % REPEAT_LENGTH];
    }
    public double getM(int m, int s, int i, int r, int j, boolean c, int diff){
        return M[m][pool[s].pattern][i][pool[r].pattern][j][c ? 1 : 0][diff + 2*minDiff];
    }
    public void setM(int m, int s, int i, int r, int j, boolean c, int diff, double val){
        //System.out.printf("%s --- %s", diff, minDiff);
        M[m][pool[s].pattern][i][pool[r].pattern][j][c ? 1 : 0][diff + 2*minDiff] = val;
    }
    public int getStrandLength(int s){
        return pool[s].repeats * REPEAT_LENGTH;
    }
    public int getNumStrands(){
        return pool.length;
    }
    public String toString(int strand){
        char[] triplet = new char[REPEAT_LENGTH * pool[strand].repeats];
        for (int i = 0; i < triplet.length; i++) triplet[i] = patternArray[pool[strand].pattern][i % REPEAT_LENGTH].toChar();
        return new String(triplet);
    }
    public int getMinDiff(){
        return minDiff;
    }
    public int getMaxDiff(){
        return maxDiff;
    }
}
