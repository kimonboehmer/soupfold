import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class SecondaryStructure {
    private class Position{
        int strand;
        int pos;
        protected Position(int strand, int pos){
            this.strand = strand;
            this.pos = pos;
        }
    }
    private int numBps;

    StrandPool sp;
    private final Position[][] partner; // -1 if unpaired
    private final int[] permutation;
    public SecondaryStructure(StrandPool sp, int m){
        this.sp = sp;
        int maxStrandLength = 0;
        for (int s = 0; s < sp.getNumStrands(); s++) if (sp.getStrandLength(s) > maxStrandLength) maxStrandLength = sp.getStrandLength(s);
        partner = new Position[m][maxStrandLength];
        for (Position[] ints : partner) Arrays.fill(ints, null);
        permutation = new int[m];
        numBps = 0;
    }

    /**
     * @param s the position of the left strand
     * @param i the position of the base on the left strand
     * @param r the position of the right strand
     * @param j the position of the base on the right strand
     */
    public void setBasePair(int s, int i, int r, int j){
        /*System.out.println("----");
        System.out.println(s);
        System.out.println(i);
        System.out.println(r);
        System.out.println(j);*/
        if (partner[s][i] != null || partner[r][j] != null) throw new RuntimeException("Position already paired!");
        partner[s][i] = new Position(r, j);
        partner[r][j] = new Position(s, i);
        numBps++;
    }
    public void setStrandRank(int s, int i){
        permutation[i] = s;
    }
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
}
