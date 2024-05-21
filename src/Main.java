import java.util.LinkedList;

public class Main {
    public static void main(String[] args) {
        testPartitionFunction();
    }
    public static void testMFE(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 27, 1);
        DP dp = new DP(sp, 3, true, new MFE());
        System.out.println((int) dp.computeMFE(3));
        System.out.println(dp.backtrack());
    }
    public static void testPartitionFunction(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 15, 1);
        DP dp = new DP(sp, 3, true, new PartitionFunction(300));
        System.out.println((int) dp.computeMFE(3));
        System.out.println(dp.backtrack());
        System.out.println(dp.dpt.E());
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
            DP dp = new DP(sp, 3, true, new MFE());
            results.add((int) dp.computeMFE(5));
        }
        System.out.println(results);
    }
    public static void testHeterodimer(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG8");
        strands.add("GUU9");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        //DP dp = new DP(sp, 3, true, new MFE());
        DP dp2 = new DP(sp, 3, true, new PartitionFunction(300));
        System.out.println(dp2.computeMFE(3));
        //System.out.println(dp.computeMFE(3));
        SecondaryStructure st = dp2.backtrack();
        System.out.println(st.toString());
        //st.toFile("heterodimer");
    }
    public static void testHeterodimer2(){
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG40");
        strands.add("GUU45");
        strands.add("CAG10");
        strands.add("GUU15");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp2 = new DP(sp, 3, true, new PartitionFunction(300));
        System.out.println(dp2.computeMFE(3));
        SecondaryStructure st = dp2.backtrack();
        System.out.println(st.toString());
    }
}
