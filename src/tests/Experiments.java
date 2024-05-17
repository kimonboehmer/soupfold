package tests;

import algorithms.DP;
import algorithms.MFE;
import algorithms.PartitionFunction;
import datastructures.Base;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;
import datastructures.TripletPool;

import java.util.LinkedList;
import java.util.Random;
public class Experiments {
    public static void testMFE(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 6, 4);
        DP dp = new DP(sp, 3, 3, true, 18, new MFE());
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
            StrandPool sp = new TripletPool(b, 27, 0);
            DP dp = new DP(sp, 5, 3, true, 81, new MFE());
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
        DP dp = new DP(sp, 3, 3, true, 30, new MFE());
        DP dp2 = new DP(sp, 3, 3, true, 30, new PartitionFunction(300));
        System.out.println(dp2.computeMFE());
        System.out.println(dp.computeMFE());
        SecondaryStructure st = dp.backtrack();
        System.out.println(st.toString());
        st.toFile("heterodimer");
    }
    public static double connectednessProbability(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 47, 2);
        DP dp = new DP(sp, 4, 3, true, 47 * 3, new PartitionFunction(300));
        DP dp2 = new DP(sp, 4, 3, false, 47 * 3, new PartitionFunction(300));
        return dp.computeMFE() / dp2.computeMFE();
    }
    public static int randomLengthTriplets(String p, int avg, int num){
        LinkedList<String> strands = new LinkedList<>();
        Random r = new Random();
        int missing = avg * num;
        for (int i = 0; i < num - 1; i++){
            int rand = r.nextInt(avg/2) + (avg-avg/4);
            missing -= rand;
            if (missing < 0) {
                strands.add(p + (rand + missing));
                break;
            }
            strands.add(p + rand);
        }
        strands.add(p + missing);
        StrandPool sp = new TripletPool(strands);
        DP dp = new DP(sp, 4, 3, true, avg*3, new MFE());
        int a = (int) dp.computeMFE();
        SecondaryStructure st = dp.backtrack();
        System.out.println(st.toString());
        return a;
    }
}
