package datastructures;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class SecondaryStructure {
    private static class Position{
        int strand;
        int pos;

        /**
         * @param strand position of the strand in the permutation
         * @param pos index in the strand
         */
        protected Position(int strand, int pos){
            this.strand = strand;
            this.pos = pos;
        }
    }

    StrandPool sp;
    private final Position[][] partner; // -1 if unpaired
    private final int[] permutation;
    private int numBPs;

    /**
     * @param sp strand pool on which the secondary structure is built
     * @param m number of interacting strands
     */
    public SecondaryStructure(StrandPool sp, int m){
        this.sp = sp;
        int maxStrandLength = 0;
        for (int s = 0; s < sp.getNumStrands(); s++) if (sp.getStrandLength(s) > maxStrandLength) maxStrandLength = sp.getStrandLength(s);
        partner = new Position[m][maxStrandLength];
        for (Position[] ints : partner) Arrays.fill(ints, null);
        permutation = new int[m];
        numBPs = 0;
    }

    /**
     * @param s the position of the left strand
     * @param i the position of the base on the left strand
     * @param r the position of the right strand
     * @param j the position of the base on the right strand
     */
    public void setBasePair(int s, int i, int r, int j){
        if (partner[s][i] != null || partner[r][j] != null) {
            throw new RuntimeException("Position already paired!");
        }
        partner[s][i] = new Position(r, j);
        partner[r][j] = new Position(s, i);
        numBPs++;
    }

    /**
     * @param s strand number
     * @param i strand position
     */
    public void setStrandRank(int s, int i){
        permutation[i] = s;
    }

    /**
     * @param s given strand
     * @param i given index
     * @return strand paired to s_i, or -1 if unpaired
     */
    public int getPairedStrand(int s, int i){
        Position p = partner[s][i];
        return p == null ? -1 : p.strand;
    }

    /**
     * @param s given strand
     * @param i given index
     * @return index (position) paired to s_i, or -1 if unpaired
     */
    public int getPairedIndex(int s, int i){
        Position p = partner[s][i];
        return p == null ? -1 : p.pos;
    }

    /**
     * @param i strand position
     * @return strand number
     */
    public int getStrandFromPosition(int i){
        return permutation[i];
    }

    /**
     * @return number of interacting strands
     */
    public int getM(){ return permutation.length; }

    /**
     * @return strand pool.
     */
    public StrandPool getStrandPool(){ return sp; }

    /**
     * @return a string representation of this secondary structure.
     */
    public String toString(){
        StringBuilder sb = new StringBuilder(sp.toString(permutation[0]));
        for (int i = 1; i < permutation.length; i++) {
            sb.append("&");
            sb.append(sp.toString(permutation[i]));
        }
        sb.append("\n");
        for (int s = 0; s < permutation.length; s++){
            for (int i = 0; i < sp.getStrandLength(permutation[s]); i++){
                Position p = partner[s][i];
                if (p == null) sb.append(".");
                else if (p.strand > s || (p.strand == s && p.pos > i)) sb.append("(");
                else sb.append(")");
            }
        }
        return sb.toString();
    }

    /**
     * @param name name of the structure, for the filename
     * stores the secondary structure in a file,
     * readable for common RNA visualization programs like "RiboSketch".
     */
    public void toFile(String name){
        try {
            File myObj = new File("structure_" + name +  ".txt");
            if (myObj.createNewFile()) {
                System.out.println("File created: " + myObj.getName());
            } else {
                System.out.println("File already exists. Overwriting...");
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        try {
            FileWriter myWriter = new FileWriter("structure_" + name +  ".sest");
            myWriter.write(toString());
            myWriter.close();
            System.out.println("Successfully wrote to the file.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    /**
     * @return the total number of base pairs in this secondary structure.
     */
    public int getNumBPs() {
        return numBPs;
    }
}
