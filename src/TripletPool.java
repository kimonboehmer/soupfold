import java.util.*;

public class TripletPool implements StrandPool {
    private static final byte REPEAT_LENGTH = 3;
    protected class Triplet {
        int pattern;
        int repeats;
        protected Triplet(int pattern, int repeats){
            this.pattern = pattern;
            this.repeats = repeats;
        }
    }
    int maxRepeats;
    int numPatterns;
    private Base[][] patternArray;
    private final Triplet[] pool;
    int[][][][][] M;
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
    }
    public TripletPool(List<String> strands){
        pool = new Triplet[strands.size()];
        maxRepeats = 0;
        HashMap<Base[], Integer> hm = new HashMap<>();
        int counter = 0;
        int i = 0;
        for (String strand : strands){
            Base[] pattern = new Base[REPEAT_LENGTH];
            for (int j = 0; j < REPEAT_LENGTH; j++){
                switch (strand.toCharArray()[j]){
                    case 'A': pattern[j] = Base.A;break;
                    case 'C': pattern[j] = Base.C;break;
                    case 'G': pattern[j] = Base.G;break;
                    case 'U': pattern[j] = Base.U;break;
                    default: throw new InputMismatchException("Invalid strand");
                }
            }
            Integer val = hm.get(pattern);
            if (val == null){
                hm.put(pattern, counter);
                val = counter++;
            }
            int reps = Integer.parseInt(strand.substring(REPEAT_LENGTH));
            if (maxRepeats < reps) maxRepeats = reps;
            Triplet t = new Triplet(val, reps);
            pool[i++] = t;
        }
        numPatterns = counter;
        patternArray = new Base[numPatterns][REPEAT_LENGTH];
        for (Map.Entry<Base[], Integer> entry : hm.entrySet()){
            patternArray[entry.getValue()] = entry.getKey();
        }
    }
    public void initializeTable(int m, int theta){
        M = new int[m+1][numPatterns][maxRepeats * REPEAT_LENGTH][numPatterns][maxRepeats * 3];
        for(int si=0;si<numPatterns;si++)for(int ii=0;ii<maxRepeats * REPEAT_LENGTH;ii++)for(int ri=0;ri<numPatterns;ri++)for(int ji=0;ji<maxRepeats * REPEAT_LENGTH;ji++){
            for (int mi = 1; mi < m+1;mi++)M[mi][si][ii][ri][ji]=DP.NOT_SET;
            M[0][si][ii][ri][ji] = 0;
        }
        for (int i = 0; i < maxRepeats * REPEAT_LENGTH; i++){
            for (int t = 0; t < Math.min(theta, maxRepeats * REPEAT_LENGTH - i); t++){
                for (int p = 0; p < numPatterns; p++){
                    M[1][p][i][p][i + t] = 0;
                }
            }
        }
    }
    public Base getBase(int s, int pos) {
        if (pos > pool[s].repeats * REPEAT_LENGTH) throw new InputMismatchException("Position greater than strand length!");
        return patternArray[pool[s].pattern][pos % REPEAT_LENGTH];
    }
    public int getM(int m, int s, int i, int r, int j){
        return M[m][pool[s].pattern][i][pool[r].pattern][j];
    }
    public void setM(int m, int s, int i, int r, int j, int val){
        M[m][pool[s].pattern][i][pool[r].pattern][j] = val;
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
