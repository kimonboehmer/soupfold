import java.util.LinkedList;

public class Main {
    public static void main(String[] args) {
        //StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 47, 1);
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG8");
        strands.add("GUU9");
        strands.add("ACG12");
        StrandPool sp = new TripletPool(strands);
        DP dp = new DP(sp, 3, true, new MFE());
        DP dp2 = new DP(sp, 3, true, new PartitionFunction(300));
        System.out.println(dp2.computeMFE(3));
        System.out.println(dp.computeMFE(3));
        SecondaryStructure st = dp.backtrack();
        System.out.println(st.toString());
        st.toFile("CAG_conn");
    }
}
