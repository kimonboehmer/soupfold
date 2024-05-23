import algorithms.DP;
import algorithms.MFE;
import algorithms.PartitionFunction;
import datastructures.Base;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;
import datastructures.TripletPool;
import experiments.Experiments;

import java.util.LinkedList;

public class Main {
    public static void main(String[] args) {
        System.out.println(connectivity());
    }
    public static void bpProbas(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 27, 1);
        for (int m = 1; m < 10; m++) System.out.println(Experiments.interactionProbability(sp, m, 100000));
    }
    public static void bpTypes(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 27, 1);
        Experiments.expNumOccurencesOfStrands(sp, 3, 10000);
    }
    public static double connectivity(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 15, 1);
        Experiments.connectivityExperiment(sp, 3, 200);
        for (int m = 1; m < 20; m++)
            System.out.printf("Not-connectedness probability with m=%d: %f\n", m, Experiments.connectivityExperiment(sp, m, 10000).getFirst());
        return 0;
    }
    public static void testMFE(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 27, 1);
        DP dp = new DP(sp, 3, 3, true, new MFE());
        System.out.println((int) dp.computeMFE());
        System.out.println(dp.backtrack());
    }
    public static void testPartitionFunction(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 27, 1);
        DP dp = new DP(sp, 3, 3, true, new PartitionFunction(300));
        System.out.println((int) dp.computeMFE());
        System.out.println(dp.backtrack());
    }
    public static void testDifferentTriplets(){
        LinkedList<Integer> results = new LinkedList<>();
        LinkedList<Base[]> ll = new LinkedList<>();
        ll.add(new Base[]{Base.C, Base.A, Base.G});
        ll.add(new Base[]{Base.A, Base.C, Base.G});
        ll.add(new Base[]{Base.G, Base.C, Base.G});
        ll.add(new Base[]{Base.G, Base.C, Base.C});
        for (Base[] b : ll){
            StrandPool sp = new TripletPool(b, 27, 1);
            DP dp = new DP(sp, 5, 3, true, new MFE());
            results.add((int) dp.computeMFE());
        }
        System.out.println(results);
    }
    public static void testHeterodimer(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG8");
        strands.add("GUU9");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp2 = new DP(sp, 3, 3, true, new PartitionFunction(300));
        System.out.println(dp2.computeMFE());
        SecondaryStructure st = dp2.backtrack();
        System.out.println(st.toString());
    }
    public static void testHeterodimer2(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG40");
        strands.add("GUU45");
        strands.add("CAG10");
        strands.add("GUU15");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp2 = new DP(sp, 3, 3, true, new PartitionFunction(300));
        System.out.println(dp2.computeMFE());
        for (int i = 0; i < 10; i++) System.out.println(dp2.backtrack());
    }
}
