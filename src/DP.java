public class DP {
    private final StrandPool sp;
    private final int INFTY = Integer.MAX_VALUE / 4;
    public static final int NOT_SET = 100000;
    double mfeValue;
    DPType dpt;
    private int startM;
    double partFuncValue;
    int theta;
    private final boolean conn;
    SecondaryStructure mfeStructure;
    public DP(StrandPool sp, int theta, boolean conn, DPType dpt){
        this.sp = sp;
        this.theta = theta;
        mfeValue = NOT_SET;
        partFuncValue = 0;
        startM = 0;
        this.conn = conn;
        this.dpt = dpt;
    }
    private double minOverStrands(int m, int r, int j){
        double res = 0;
        for (int t = 0; t < sp.getNumStrands();t++){
            double val = getOrComputeM(m, t, 0, r, j, conn);
            res = dpt.min(res, val);
        }
        return res;
    }
    private double minOverSecStrands(int m, int s, int i){
        double res = 0;
        for (int u = 0; u < sp.getNumStrands();u++){
            double val = getOrComputeM(m, s, i, u, sp.getStrandLength(u) - 1, conn);
            res = dpt.min(res, val);
        }
        return res;
    }
    private double M(int m, int s, int i, int r, int j, boolean c){
        if (m==0) return dpt.noEffect();
        if (m==1) {
            if (i+1 < sp.getStrandLength(s) && j > 0 && j > i+1) return getOrComputeM(1, s, i+1, s, j-1, false);
            return 0;
        }
        if (j > 0){
            if (i+1 < sp.getStrandLength(s)) return getOrComputeM(m, s, i+1, r, j-1, c);
            if (m==2) return getOrComputeM(1, r, 0, r, j-1, false);
            return minOverStrands(m-1, r, j-1);
        }
        if (i+1 < sp.getStrandLength(s)){
            if (m==2) return getOrComputeM(1, s, i+1, s, sp.getStrandLength(s)-1, false);
            return minOverSecStrands(m-1, s, i+1);
        }
        return 0;
    }
    private boolean bp(int s, int i, int r, int j){
        return (Base.pair(sp.getBase(s, i), sp.getBase(r, j)));
    }
    private double computeMultiloop(int m, int s, int i, int r, int j, boolean c){
        double mfe = 0;
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s);
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k < lim; k++) if (bp(s,i,r,j)){
                    double val = dpt.sum(dpt.sum(dpt.E(), M(1, s, i, s, k, false)), M(m, s, k, r, j+1, (m > 1) && c));
                    mfe = dpt.min(mfe, val);
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++) if (bp(s,i,r,j)){
                    double val = dpt.sum(dpt.sum(dpt.E(), M(m, s, i, r, k, false)),M(1, r, k, r, j+1, false));
                    mfe = dpt.min(mfe, val);
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                for (int k = 0; k < sp.getStrandLength(t); k++) if(bp(s, i, t, k)){
                    double val = dpt.sum(dpt.sum(dpt.E(), M(mm, s, i, t, k, false)), M(m-mm+1, t, k, r, j+1, c));
                    mfe = dpt.min(mfe, val);
                }
            }
        }
        return mfe;
    }

    /**
     * @param m number of still available strands
     * @param s index of leftmost strand
     * @param i leftmost position on s
     * @param r index of rightmost strand
     * @param j rightmost position on r
     * @param c true if s and r have to be connected, false else
     * @return minimum free energy value for an interval from s_i to r_j on m strands, respecting the connectivity bit c.
     */
    private double computeMFERegion(int m, int s, int i, int r, int j, boolean c){
        double unpaired;
        if (i+1 < sp.getStrandLength(s)) unpaired = getOrComputeM(m, s, i+1, r, j, c);
        else if (m==2) unpaired = getOrComputeM(1, r, 0, r, j, false);
        else unpaired = minOverStrands(m - 1, r, j);
        double multi = computeMultiloop(m, s, i, r, j, c);
        double mfe = dpt.min(unpaired, multi);
        sp.setM(m, s, i, r, j, c, mfe);
        return mfe;
    }
    private double getOrComputeM(int m, int s, int i, int r, int j, boolean c){
        double val = sp.getM(m, s, i, r, j, c);
        if (val == NOT_SET){
            val = computeMFERegion(m, s, i, r, j, c);
        }
        return val;
    }
    public double computeMFE(int m){
        startM = m;
        sp.initializeTable(m, theta, dpt.initValue());
        double mfe = 0;
        for (int s = 0; s < sp.getNumStrands(); s++) {
            for (int r = 0; r < sp.getNumStrands(); r++) {
                double val = getOrComputeM(m, s, 0, r, sp.getStrandLength(r) - 1, conn);
                mfe = dpt.min(mfe, val);
            }
        }
        mfeValue = mfe;
        return mfe;
    }

    /**
     * assumes that the DP table is already filled (i.e. computeMFE already executed)!
     * @return one optimal secondary structure
     */
    public SecondaryStructure backtrack(){
        mfeStructure = new SecondaryStructure(sp, startM);
        for (int s = 0; s < sp.getNumStrands(); s++){
        for (int r = 0; r < sp.getNumStrands(); r++) {
            if (sp.getM(startM, s, 0, r, sp.getStrandLength(r) - 1, conn) == mfeValue) {
                mfeStructure.setStrandRank(s, 0);
                mfeStructure.setStrandRank(r, startM - 1);
                recBacktrack(s, 0, r, sp.getStrandLength(r) - 1,  conn, 0, startM - 1, mfeValue);
                return mfeStructure;
            }
        }
        }
        return mfeStructure;
    }
    private boolean recBacktrack(int s, int i, int r, int j, boolean c, int left, int right, double val) {
        System.out.println("-------------");
        System.out.printf("New Backtrack with: %s, %s, %s, %s", left, i, right, j);
        System.out.println("-------------");
        if (val == 0) return true;
        if (left == right && j - i <= theta) return true;
        int m = right - left + 1;
        // unpaired
        if (i + 1 < sp.getStrandLength(s)){
            if (sp.getM(m, s, i + 1, r, j, c) == val) return recBacktrack(s, i + 1, r, j, c, left, right, val);
        }
        else{
            if (m == 2) {
                if (sp.getM(1, r, 0, r, j, false) == val) return recBacktrack(r, 0, r, j, false, right, right, val);
            }
            else for (int t = 0; t < sp.getNumStrands(); t++)
            if (sp.getM(m - 1, t, 0, r, j, conn) == val) {
                mfeStructure.setStrandRank(t, left+1);
                return recBacktrack(t, 0, r, j, conn,left + 1, right, val);
            }
        }
        // multiloop
        return backtrackMultiloop(m, s, i, r, j, c, left, right, val);
    }
    private boolean backtrackMultiloop(int m, int s, int i, int r, int j, boolean c, int left, int right, double val){
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s);
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k < lim; k++){
                    double m1 = M(1, s, i, s, k, false);
                    double m2 = M(m, s, k, r, j+1, m>1 && c);
                    if (bp(s,i,s,k) && val == dpt.sum(dpt.sum(dpt.E(), m1), m2)){
                        mfeStructure.setBasePair(left, i, left, k);
                        return     backtrackHelper(m, s, i, s, k, false, left, left, m1)
                                && backtrackHelper(m, s, k, r, j+1, (m>1 && c), left, right, m2);
                }
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++){
                    double m1 = M(m, s, i, r, k, false);
                    double m2 = M(1, r, k, r, j+1, false);
                    if (bp(s,i,r,k) && val == dpt.sum(dpt.sum(dpt.E(), m1),m2)){
                        mfeStructure.setBasePair(left, i, right, k);
                        return     backtrackHelper(m, s, i, r, k, false, left, right, m1)
                                && backtrackHelper(1, r, k, r, j+1, false, right, right, m2);
                }}
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                    for (int k = 0; k < sp.getStrandLength(t); k++){
                        double m1 = M(mm, s, i, t, k, false);
                        double m2 = M(m-mm+1, t, k, r, j+1, c);
                        if (bp(s,i,t,k) && val == dpt.sum(dpt.sum(dpt.E(),m1), m2)){
                            mfeStructure.setBasePair(left, i, left+mm-1, k);
                            mfeStructure.setStrandRank(t, left + mm - 1);
                            return     backtrackHelper(mm, s, i, t, k, false, left, left + mm - 1, m1)
                                    && backtrackHelper(m-mm+1, t, k, r, j+1, c, left + mm - 1, right, m2);
                    }}
                }
        }
        System.out.println("Cannot find the follow-up DP decision");
        return false;
    }
    private boolean backtrackHelper(int m, int s, int i, int r, int j, boolean c, int left, int right, double val){
        if (m == 1){
            if (j - i < 3) return true; // empty region
            if (j > 0 && i+1 < sp.getStrandLength(s) && sp.getM(1, s, i+1, s, j-1, false) == val){
            return recBacktrack(s, i+1, s, j-1, false, left, right, val);
        }
        }
        else if (j > 0){
            if (i+1 < sp.getStrandLength(s)){
                if (sp.getM(m, s, i+1, r, j-1, c) == val){
                    return recBacktrack(s, i+1, r, j-1, c, left, right, val);
                }
            }
            else if (m == 2){
                if (sp.getM(1, r, 0, r, j-1, false) == val){
                    return recBacktrack(r, 0, r, j-1, false, right, right, sp.getM(1, r, 0, r, j-1, false));
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                double m1 = sp.getM(m-1, t, 0, r, j-1, conn);
                if (m1 == val){
                        mfeStructure.setStrandRank(t, left+1);
                        return recBacktrack(t, 0, r, j-1, conn,left+1, right, m1);
                    }
                }
        }
        else{
            if (i+1 < sp.getStrandLength(s)){
                if (m==2) {
                    if (sp.getM(1, s, i + 1, s, sp.getStrandLength(s) - 1, false) == val) {
                        return recBacktrack(s, i + 1, s, sp.getStrandLength(s) - 1, false, left, left, sp.getM(1, s, i + 1, s, sp.getStrandLength(s) - 1, false));
                    }
                }
                else for (int u = 0; u < sp.getNumStrands(); u++){
                    double m1 = sp.getM(m-1, s, i, u, sp.getStrandLength(u)-1, conn);
                    if (m1 == val){
                        mfeStructure.setStrandRank(u, right-1);
                        return recBacktrack(s, i, u, sp.getStrandLength(u)-1, conn, left, right-1, m1);
                    }
                }
            }
            return m <= 2; //empty region
        }
        //System.out.printf("!!!!!!! -- %s %s %s %s -- ", left, i, right, j);
        return false;
    }
}
