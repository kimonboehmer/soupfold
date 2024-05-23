package datastructures;

import algorithms.DP;
import datastructures.Base;
import datastructures.StrandPool;

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
    double[][][][][][] M;
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
            char[] strandCharArray = strand.toCharArray();
            Base[] pattern = new Base[REPEAT_LENGTH];
            for (int j = 0; j < REPEAT_LENGTH; j++){
                pattern[j] = Base.toBase(strandCharArray[j]);
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
    public void initializeTable(int m, int theta, double initValue){
        M = new double[m+1][numPatterns][maxRepeats * REPEAT_LENGTH][numPatterns][maxRepeats * 3][2];
        for(int si=0;si<numPatterns;si++)for(int ii=0;ii<maxRepeats * REPEAT_LENGTH;ii++)for(int ri=0;ri<numPatterns;ri++)for(int ji=0;ji<maxRepeats * REPEAT_LENGTH;ji++){
            for (int mi = 1; mi < m+1;mi++){
                M[mi][si][ii][ri][ji][0]= DP.NOT_SET;
                M[mi][si][ii][ri][ji][1]= DP.NOT_SET;
            }
            M[0][si][ii][ri][ji][0] = initValue;
        }
        for (int i = 0; i < REPEAT_LENGTH; i++){
            for (int t = 0; t <= theta; t++){
                for (int p = 0; p < numPatterns; p++){
                    M[1][p][i][p][t][0] = initValue;
                }
            }
        }
    }
    public Base getBase(int s, int pos) {
        if (pos > pool[s].repeats * REPEAT_LENGTH) throw new InputMismatchException("Position greater than strand length!");
        return patternArray[pool[s].pattern][pos % REPEAT_LENGTH];
    }
    public double getM(int m, int s, int i, int r, int j, boolean c){
        if (m == 1) return M[1][pool[s].pattern][i % REPEAT_LENGTH][pool[s].pattern][j - i][0];
        return M[m][pool[s].pattern][getStrandLength(s) - i - 1][pool[r].pattern][j][c ? 1 : 0];
    }
    public void setM(int m, int s, int i, int r, int j, boolean c, double val){
        if (m == 1) M[1][pool[s].pattern][i % REPEAT_LENGTH][pool[s].pattern][j - i][0] = val;
        M[m][pool[s].pattern][getStrandLength(s) - i - 1][pool[r].pattern][j][c ? 1 : 0] = val;
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
}
