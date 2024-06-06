package datastructures;
public interface StrandPool {
    /**
     * @return the total number of different strands in the soup.
     */
int getNumStrands();
    /**
     * @param m number of remaining strands
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if structure is connected
     * @return MFE/PF value that is stored in the DP table
     */
double getM(int m, int s, int i, int r, int j, boolean c);

    /**
     * @param m number of remaining strands
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if structure is connected
     * @param value MFE/PF value that should be set in the DP table
     */
void setM(int m, int s, int i, int r, int j, boolean c, double value);

    /**
     * @param s strand number
     * @param pos position on s
     * @return a Base (A,C,G or U)
     */
Base getBase(int s, int pos);

    /**
     * @param s strand number
     * @return number of bases in s
     */
int getStrandLength(int s);
    /**
     * @param m number of interacting strands
     * @param theta minimum base pair span
     * @param initValue value for an empty structure
     */
void initializeTable(int m, int theta, double initValue);

    /**
     * @param strand strand number
     * @return its representation as a string
     */
    default String toString(int strand) {
        char[] c = new char[getStrandLength(strand)];
        for (int i = 0; i < getStrandLength(strand); i++) c[i] = getBase(strand, i).toChar();
        return new String(c);
    }

}
