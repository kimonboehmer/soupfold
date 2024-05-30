package experiments;

import datastructures.Base;

import java.io.*;
import java.util.ArrayList;
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
    public static double[][] changeZeros(double[][] data, double[][] otherData){
        for (int i = 0; i < 64; i++) data[i][i] = otherData[i][i];
        return data;
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
    public static void createGraph(double[][] data, String path) throws IOException {
        StringBuilder sb = new StringBuilder("#64 100");
        for (int i = 0; i < 64; i++) for (int j = i+1; j < 64; j++){
            if (data[i][i] < 0 || data[j][j] < 0 || data[i][j] < 0.35){
                sb.append("\n").append(i).append(" ").append(j);
            }
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(path));
        writer.write(sb.toString());
        writer.close();
    }
    public static void getClique(){
        boolean[] out = new boolean[64];
        int[] arr = new int[]{42,40,34,32,31,29,23,21,20,63,61,8,5,4,2,1,0,17,16,55,10,53,43,62,59,58,47,46,27,14,57,44,39,35,30,56,11,54,50,45,41,38,36,26,60,6,19,15,51,49,33,28,24,9,7,18,13,37,25,3,12};
        for (int el : arr) out[el] = true;
        int i = 0;
        for (Base X : bases) for (Base Y : bases) for (Base Z : bases){
            if (!out[i++]) System.out.println(tripletToString(new Base[]{X,Y,Z}));
        }
    }
    public static String tripletToString(Base[] triplet){
        return new String(new char[]{triplet[0].toChar(), triplet[1].toChar(), triplet[2].toChar()});
    }
}

