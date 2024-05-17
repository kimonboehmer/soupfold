import java.util.LinkedList;

public class Main {
    public static void main(String[] args) {
        testMFE();
    }
    public static void testMFE(){
        StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 6, 4);
        DP dp = new DP(sp, 3, true, 18, new MFE());
        System.out.println((int) dp.computeMFE(2));
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
            DP dp = new DP(sp, 3, true, 81, new MFE());
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
        DP dp = new DP(sp, 3, true, 30, new MFE());
        DP dp2 = new DP(sp, 3, true, 30, new PartitionFunction(300));
        System.out.println(dp2.computeMFE(3));
        System.out.println(dp.computeMFE(3));
        SecondaryStructure st = dp.backtrack();
        System.out.println(st.toString());
        st.toFile("heterodimer.secstr");
    }
}
