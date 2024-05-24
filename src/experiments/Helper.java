package experiments;

import datastructures.Base;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

import static datastructures.Base.*;
import static datastructures.Base.U;
import static experiments.FinalExperiments.bases;

public class Helper {
    public static double[][] readCSV(String path){
        List<double[]> records = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] values = line.split(";");
                double[] vals = new double[values.length - 1];
                for (int i = 1; i < values.length; i++) {
                    if (Objects.equals(values[i], "")) vals[i-1] = 0;
                    else vals[i-1] = Double.parseDouble(values[i]);
                }
                records.add(vals);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return records.toArray(double[][]::new);
    }
    public static double[][] changeZeros(double[][] data){
        double[][] newData = new double[64][64];
        int i = -1;
        for (Base T : bases) for (Base V : bases) for (Base W : bases) {
            i++;
            int j = -1;
            int diff = 0;
            for (Base X : bases) for (Base Y : bases) for (Base Z : bases) {
                j++;
                HashSet<Base> h = new HashSet<>(List.of(new Base[]{T, V, W, X, Y, Z}));
                if ( (X.ordinal() < T.ordinal() || (X.ordinal() == T.ordinal() && Y.ordinal() < V.ordinal()) || (X.ordinal() == T.ordinal() && Y.ordinal() == V.ordinal() && Z.ordinal() <= W.ordinal()))){
                    newData[i][j] = -2;
                    diff++;
                }
                else if (h.size() == 1 || (h.size() == 2 && !Base.pair((Base) h.toArray()[0], (Base) h.toArray()[1]))){
                    newData[i][j] = -1;
                    diff++;
                }
                else{
                    newData[i][j] = data[i][j - diff];
                }
            }
        }
        return newData;
    }
    public static void writeCSV(double[][] results, String path) throws IOException {
        StringBuilder sb = new StringBuilder(";");
        String[] triplets = tripletsInOrder();
        for (int i = 0; i < 64; i++){
            sb.append(triplets[i]);
            if(i < 63) sb.append(";");
        }
        sb.append("\n");
        for(int i = 0; i < 64; i++) {
            sb.append(triplets[i]);
            for(int j = 0; j < 64; j++) {
                sb.append(";");
                sb.append(results[i][j]);
            }
            sb.append("\n");
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(path));
        writer.write(sb.toString());
        writer.close();
    }
    public static String[] tripletsInOrder(){
        String[] s = new String[64];
        Base[] bases = new Base[]{A, C, G, U};
        int i = 0;
        for (Base X : bases) for (Base Y : bases) for (Base Z : bases){
            s[i++] = new String(new char[]{X.toChar(), Y.toChar(), Z.toChar()});
        }
        return s;
    }
}
