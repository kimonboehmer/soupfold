import java.util.LinkedList;

public class Main {
    public static void main(String[] args) {
        //StrandPool sp = new TripletPool(new Base[]{Base.C, Base.A, Base.G}, 47, 1);
        LinkedList<String> strands = new LinkedList<>();
        strands.add("CAG47");
        strands.add("GUU50");
        strands.add("ACG45");
        StrandPool sp = new TripletPool(strands);
        DP dp = new DP(sp, 3);
        System.out.println(dp.computeMFE(3));
        SecondaryStructure st = dp.backtrack();
        System.out.println(st.toString());
        st.toFile("CAG");
    }
}
