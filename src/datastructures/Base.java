package datastructures;

import java.util.InputMismatchException;

public enum Base {
    /**
     * adenine
     */
    A,
    /**
     * cytosine
     */
    C,
    /**
     * guanine
     */
    G,
    /**
     * uracil
     */
    U;
    private static final boolean[][] BASEPAIRS = new boolean[][]{
            new boolean[]{false, false, false, true},
            new boolean[]{false, false, true, false},
            new boolean[]{false, true, false, true},
            new boolean[]{true, false, true, false}};
    private static final char[] CHARS = new char[]{'A','C','G','U'};

    /**
     * @param x first letter
     * @param y second letter
     * @return true if the two letters are Watson-Crick or Wobble pairs
     */
    public static boolean pair(Base x, Base y){
        return BASEPAIRS[x.ordinal()][y.ordinal()];
    }

    /**
     * @return the character associated to this base.
     */
    public char toChar(){
        return CHARS[this.ordinal()];
    }

    /**
     * @param c a character from A,C,G,U
     * @return the corresponding base
     */
    public static Base toBase(char c){
        return switch (c) {
            case 'A' -> Base.A;
            case 'C' -> Base.C;
            case 'G' -> Base.G;
            case 'U' -> Base.U;
            default -> throw new InputMismatchException("Invalid strand");
        };
    }
}
