package experiments;
import algorithms.DP;
import algorithms.MFE;
import algorithms.PartitionFunction;
import datastructures.Base;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;
import datastructures.TripletPool;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static datastructures.Base.*;
import static experiments.Sampling.avgBPs;
import static experiments.Sampling.classifyBasePairs;

public class FinalExperiments {
    public static Base[] bases = new Base[]{A, C, G, U};
    public static void bpTypes() throws IOException {
        double[][][][] results = new double[64][64][4][3];
        int cV = 0;
        for (Base T : bases) for (Base V : bases) for (Base W : bases) {
            int cH = 0;
            for (Base X : bases) for (Base Y : bases) for (Base Z : bases) {
                if (X.ordinal() < T.ordinal() || (X.ordinal() == T.ordinal() && Y.ordinal() < V.ordinal()) || (X.ordinal() == T.ordinal() && Y.ordinal() == V.ordinal() && Z.ordinal() <= W.ordinal())){
                    cH++;
                    continue;
                }
                String triplet1 = new String(new char[]{T.toChar(), V.toChar(), W.toChar()});
                String triplet2 = new String(new char[]{X.toChar(), Y.toChar(), Z.toChar()});
                LinkedList<String> st = new LinkedList<>();
                st.add(triplet1 + "20");
                st.add(triplet2 + "20");
                System.out.println("Triplet " + triplet1 + " and " + triplet2);
                HashSet<Base> h = new HashSet<>(List.of(new Base[]{T, V, W, X, Y, Z}));
                if (h.size() == 1 || (h.size() == 2 && !Base.pair((Base) h.toArray()[0], (Base) h.toArray()[1]))) {
                    System.out.println("No base pairs possible");
                    for (int m = 0; m < 4; m++) for (int t = 0; t < 3; t++) results[cV][cH][m][t] = -1;
                    cH++;
                    continue;
                }
                for (int m = 2; m < 6; m++) {
                    TripletPool tp = new TripletPool(st);
                    double[] data = classifyBasePairs(tp, m, 1000);
                    System.out.printf("For m=%d, interior bps: %f, exterior homo-bps: %f, exterior hetero-bps: %f\n", m, data[0], data[1], data[2]);
                    System.arraycopy(data, 0, results[cV][cH][m - 2], 0, 3);
                }
                cH++;
            }
            cV++;
        }
        for (int m = 0; m < 4; m++) for (int t = 0; t < 3; t++){
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < 64; i++) {
                for(int j = 0; j < 64; j++) {
                    sb.append(results[i][j][m][t]);
                    if(j < 63)
                        sb.append(";");
                }
                sb.append("\n");
            }
            BufferedWriter writer = new BufferedWriter(new FileWriter("table_m"+(m+2)+"_t"+t+".csv"));
            writer.write(sb.toString());
            writer.close();
        }
    }
    public static void bpProbas(){
        StrandPool sp = new TripletPool(new Base[]{C, Base.A, Base.G}, 27, 1);
        for (int m = 1; m < 10; m++) System.out.println(Sampling.interactionProbability(sp, m, 100000));
    }
    public static void strandTypes(){
        StrandPool sp = new TripletPool(new Base[]{C, Base.A, Base.G}, 27, 1);
        Sampling.expNumOccurencesOfStrands(sp, 3, 10000);
    }
    public static double connectivity(){
        StrandPool sp = new TripletPool(new Base[]{C, Base.A, Base.G}, 15, 1);
        for (int m = 1; m < 20; m++)
            System.out.printf("Not-connectedness probability with m=%d: %f\n", m, Sampling.connectivityExperiment(sp, m, 10000).getFirst());
        return 0;
    }
    public static void testMFE(){
        for (Base X : bases) for (Base Y : bases) for (Base Z : bases) {
            StrandPool sp = new TripletPool(new Base[]{X,Y,Z}, 20, 1);
            DP dp = new DP(sp, 3, 3, false, new MFE());
            System.out.println(X.toChar());
            System.out.println(Y.toChar());
            System.out.println(Z.toChar());
            System.out.println((int) dp.getMFE());
            System.out.println("__________");
        }
    }
    public static void MFEStructure(){
        StrandPool sp = new TripletPool(new Base[]{C,A,G}, 3, 1);
        DP dp = new DP(sp, 4, 3, false, new MFE());
        System.out.println(dp.backtrack());
    }
    public static void MFEIncreasingStrandLength(){
        for (int n = 1; n < 50; n++){
            StrandPool sp = new TripletPool(new Base[]{C, Base.A, Base.G}, n, 0);
            DP dp = new DP(sp, 4, 3, true, new MFE());
            System.out.printf("n=%d: %d\n", n, (int) dp.getMFE());
        }
    }
    public static void MFEIncreasingStrandNumber(){
        for (int m = 1; m < 20; m++){
            StrandPool sp = new TripletPool(new Base[]{C, Base.A, Base.G}, 15, 0);
            DP dp = new DP(sp, m, 3, true, new MFE());
            System.out.printf("m=%d: %d\n", m, (int) dp.getMFE());
        }
    }
    public static void testPartitionFunction(){
        StrandPool sp = new TripletPool(new Base[]{Base.G, C, Base.G}, 47, 30, 5);
        for (int i = 0; i < 5 ; i++) System.out.println(sp.getStrandLength(i));
        DP dp = new DP(sp, 4, 3, true, new PartitionFunction(300));
        System.out.println((int) dp.getMFE());
        System.out.println(dp.backtrack());
    }
    public static void testDifferentTriplets(){
        LinkedList<Integer> results = new LinkedList<>();
        LinkedList<Base[]> ll = new LinkedList<>();
        ll.add(new Base[]{C, Base.A, Base.G});
        ll.add(new Base[]{Base.A, C, Base.G});
        ll.add(new Base[]{Base.G, C, Base.G});
        ll.add(new Base[]{Base.G, C, C});
        for (Base[] b : ll){
            StrandPool sp = new TripletPool(b, 47, 1);
            DP dp = new DP(sp, 4, 3, true, new MFE());
            results.add((int) dp.getMFE());
        }
        System.out.println(results);
    }
    public static void avgBPForDifferentTripletPatterns(){
        LinkedList<Double> results = new LinkedList<>();
        LinkedList<Base[]> ll = new LinkedList<>();
        LinkedList<Integer> sizes = new LinkedList<>();
        Random rand = new Random();
        for (int i = 0; i < 5; i++) sizes.add(rand.nextInt(20, 70));
        ll.add(new Base[]{C, Base.A, Base.G});
        ll.add(new Base[]{Base.A, C, Base.G});
        ll.add(new Base[]{Base.G, C, Base.G});
        ll.add(new Base[]{Base.G, C, C});
        for (Base[] b : ll){
            String s = new String(new char[]{b[0].toChar(), b[1].toChar(), b[2].toChar()});
            LinkedList<String> str = new LinkedList<>();
            for (int size : sizes) str.add(s + size);
            StrandPool sp = new TripletPool(str);
            double d = avgBPs(sp, 4, 10000);
            results.add(d);
            System.out.printf("Avg base pairs of triplet pattern %s: %f\n", s, d);
        }
        System.out.println(results);
    }
    public static void maxBPForDifferentTripletPatterns(){
        LinkedList<Double> results = new LinkedList<>();
        LinkedList<Base[]> ll = new LinkedList<>();
        LinkedList<Integer> sizes = new LinkedList<>();
        Random rand = new Random();
        for (int i = 0; i < 5; i++) sizes.add(rand.nextInt(20, 70));
        ll.add(new Base[]{C, Base.A, Base.G});
        ll.add(new Base[]{Base.A, C, Base.G});
        ll.add(new Base[]{Base.G, C, Base.G});
        ll.add(new Base[]{Base.G, C, C});
        for (Base[] b : ll){
            String s = new String(new char[]{b[0].toChar(), b[1].toChar(), b[2].toChar()});
            LinkedList<String> str = new LinkedList<>();
            for (int size : sizes) str.add(s + size);
            StrandPool sp = new TripletPool(str);
            DP dp = new DP(sp, 4, 3, true, new MFE());
            results.add(dp.getMFE());
            System.out.printf("Max base pairs of triplet pattern %s: %f\n", s, dp.getMFE());
            System.out.println(dp.backtrack());
        }
        System.out.println(results);
    }
    public static void testHeteroTriplets(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG8");
        strands.add("GUU9");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp2 = new DP(sp, 3, 3, true, new MFE());
        System.out.println(dp2.getMFE());
        SecondaryStructure st = dp2.backtrack();
        System.out.println(st.toString());
    }
    public static void testHeteroTriplets2(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG40");
        strands.add("GUU45");
        strands.add("CAG10");
        strands.add("GUU15");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp2 = new DP(sp, 3, 3, true, new PartitionFunction(300));
        System.out.println(dp2.getMFE());
        for (int i = 0; i < 10; i++) System.out.println(dp2.backtrack());
    }
    public static void doExperiment(int num) throws IOException {
        switch (num) {
            case 0 -> bpProbas();
            case 1 -> bpTypes();
            case 2 -> connectivity();
            case 3 -> testMFE();
            case 4 -> testPartitionFunction();
            case 5 -> testDifferentTriplets();
            case 6 -> testHeteroTriplets();
            case 7 -> testHeteroTriplets2();
            case 8 -> avgBPForDifferentTripletPatterns();
            case 9 -> maxBPForDifferentTripletPatterns();
            case 10 -> strandTypes();
            case 11 -> MFEIncreasingStrandLength();
            case 12 -> MFEIncreasingStrandNumber();
            case 13 -> MFEStructure();
        }
    }
}
