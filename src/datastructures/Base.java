package datastructures;

import java.util.InputMismatchException;

public enum Base {
    A,C,G,U;
    private static final boolean[][] BASEPAIRS = new boolean[][]{
            new boolean[]{false, false, false, true},
            new boolean[]{false, false, true, false},
            new boolean[]{false, true, false, true},
            new boolean[]{true, false, true, false}};
    private static final char[] CHARS = new char[]{'A','C','G','U'};
    public static boolean pair(Base x, Base y){
        return BASEPAIRS[x.ordinal()][y.ordinal()];
    }
    public char toChar(){
        return CHARS[this.ordinal()];
    }
    public static Base toBase(char c){
        switch (c){
            case 'A': return Base.A;
            case 'C': return Base.C;
            case 'G': return Base.G;
            case 'U': return Base.U;
            default: throw new InputMismatchException("Invalid strand");
        }
    }
}
