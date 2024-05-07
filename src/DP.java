public class DP {
    private final StrandPool sp;
    private final int INFTY = Integer.MAX_VALUE / 4;
    public static final int NOT_SET = 100000;
    int mfeValue;
    private int startM;
    int partFuncValue;
    int theta;
    SecondaryStructure mfeStructure;
    public DP(StrandPool sp, int theta){
        this.sp = sp;
        this.theta = theta;
        mfeValue = NOT_SET;
        partFuncValue = 0;
        startM = 0;
    }
    private int minOverStrands(int m, int r, int j){
        int res = 0;
        for (int t = 0; t < sp.getNumStrands();t++){
            int val = getOrComputeM(m, t, 0, r, j);
            if (val < res) res = val;
        }
        return res;
    }
    private int minOverSecStrands(int m, int s, int i){
        int res = 0;
        for (int u = 0; u < sp.getNumStrands();u++){
            int val = getOrComputeM(m, s, i, u, sp.getStrandLength(u) - 1);
            if (val < res) res = val;
        }
        return res;
    }
    private int M(int m, int s, int i, int r, int j){
        if (m==0) return 0;
        if (m==1) {
            if (i+1 < sp.getStrandLength(s) && j > 0 && j > i+1) return getOrComputeM(1, s, i+1, s, j-1);
            return 0;
        }
        if (j > 0){
            if (i+1 < sp.getStrandLength(s)) return getOrComputeM(m, s, i+1, r, j-1);
            if (m==2) return getOrComputeM(1, r, 0, r, j-1);
            return minOverStrands(m-1, r, j-1);
        }
        if (i+1 < sp.getStrandLength(s)){
            if (m==2) return getOrComputeM(1, s, i+1, s, sp.getStrandLength(s)-1);
            return minOverSecStrands(m-1, s, i+1);
        }
        return 0;
    }
    private int E(int s, int i, int r, int j){
        if (Base.pair(sp.getBase(s, i), sp.getBase(r, j))) return -1;
        return INFTY;
    }
    private int computeMultiloop(int m, int s, int i, int r, int j){
        int mfe = INFTY;
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s);
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k < lim; k++){
                    int val = E(s, i, s, k) + M(1, s, i, s, k) + M(m, s, k, r, j+1);
                    if (val < mfe) mfe = val;
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++){
                    int val = E(s, i, r, k) + M(m, s, i, r, k) + M(1, r, k, r, j+1);
                    if (val < mfe) mfe = val;
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                for (int k = 0; k < sp.getStrandLength(t); k++){
                    int val = E(s, i, t, k) + M(mm, s, i, t, k) + M(m-mm+1, t, k, r, j+1);
                    if (val < mfe) mfe = val;
                }
            }
        }
        return mfe;
    }
    private int computeMFERegion(int m, int s, int i, int r, int j){
        int unpaired;
        if (i+1 < sp.getStrandLength(s)) unpaired = getOrComputeM(m, s, i+1, r, j);
        else if (m==2) unpaired = getOrComputeM(1, r, 0, r, j);
        else unpaired = minOverStrands(m - 1, r, j);
        int stack = 0;//E(s, i, r, j) + M(m, s, i, r, j);
        int multi = computeMultiloop(m, s, i, r, j);
        int mfe = min(unpaired, stack, multi);
        sp.setM(m, s, i, r, j, mfe);
        return mfe;
    }
    private int getOrComputeM(int m, int s, int i, int r, int j){
        int val = sp.getM(m, s, i, r, j);
        if (val == NOT_SET){
            val = computeMFERegion(m, s, i, r, j);
        }
        return val;
    }
    public int computeMFE(int m){
        startM = m;
        sp.initializeTable(m, theta);
        int mfe = INFTY;
        for (int s = 0; s < sp.getNumStrands(); s++) {
            for (int r = 0; r < sp.getNumStrands(); r++) {
                int val = getOrComputeM(m, s, 0, r, sp.getStrandLength(r) - 1);
                if (val < mfe) mfe = val;
            }
        }
        mfeValue = mfe;
        return mfe;
    }
    private static int min(int a, int b, int c){
        if (b < a) a = b;
        if (c < a) a = c;
        return a;
    }

    /**
     * assumes that the DP table is already filled (i.e. computeMFE already executed)!
     * @return one optimal secondary structure
     */
    public SecondaryStructure backtrack(){
        mfeStructure = new SecondaryStructure(sp, startM);
        for (int s = 0; s < sp.getNumStrands(); s++){
        for (int r = 0; r < sp.getNumStrands(); r++) {
            if (sp.getM(startM, s, 0, r, sp.getStrandLength(r) - 1) == mfeValue) {
                mfeStructure.setStrandRank(s, 0);
                mfeStructure.setStrandRank(r, startM - 1);
                recBacktrack(s, 0, r, sp.getStrandLength(r) - 1, 0, startM - 1, mfeValue);
                return mfeStructure;
            }
        }
        }
        return mfeStructure;
    }
    private boolean recBacktrack(int s, int i, int r, int j, int left, int right, int val) {
        System.out.println("-------------");
        System.out.printf("New Backtrack with: %s, %s, %s, %s", left, i, right, j);
        System.out.println("-------------");
        if (val == 0) return true;
        if (left == right && j - i <= theta) return true;
        int m = right - left + 1;
        // unpaired
        if (i + 1 < sp.getStrandLength(s)){
            if (sp.getM(m, s, i + 1, r, j) == val) return recBacktrack(s, i + 1, r, j, left, right, val);
        }
        else{
            if (m == 2) {
                if (sp.getM(m - 1, r, 0, r, j) == val) return recBacktrack(r, 0, r, j, right, right, val);
            }
            else for (int t = 0; t < sp.getNumStrands(); t++)
            if (sp.getM(m - 1, t, 0, r, j) == val) {
                mfeStructure.setStrandRank(t, left+1);
                return recBacktrack(t, 0, r, j, left + 1, right, val);
            }
        }
        // stack
        /*if (backtrackHelper(m, s, i, r, j, left, right, val - E(s, i, r, j))) {
            mfeStructure.setBasePair(left, i, right, j);
            return true;
        }*/
        // multiloop
        return backtrackMultiloop(m, s, i, r, j, left, right, val);
    }
    private boolean backtrackMultiloop(int m, int s, int i, int r, int j, int left, int right, int val){
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s);
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k < lim; k++){
                    int cand = E(s, i, s, k) + M(1, s, i, s, k) + M(m, s, k, r, j+1);
                    if (val == cand){
                        mfeStructure.setBasePair(left, i, left, k);
                        return     backtrackHelper(m, s, i, s, k, left, left, M(1, s, i, s, k))
                                && backtrackHelper(m, s, k, r, j+1, left, right, M(m, s, k, r, j+1));
                }
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++){
                    int cand = E(s, i, r, k) + M(m, s, i, r, k) + M(1, r, k, r, j+1);
                    if (val == cand){
                        mfeStructure.setBasePair(left, i, right, k);
                        return     backtrackHelper(m, s, i, r, k, left, right, M(m, s, i, r, k))
                                && backtrackHelper(1, r, k, r, j+1, right, right, M(1, r, k, r, j+1));
                }}
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                    for (int k = 0; k < sp.getStrandLength(t); k++){
                        int cand = E(s, i, t, k) + M(mm, s, i, t, k) + M(m-mm+1, t, k, r, j+1);
                        if (val == cand){
                            mfeStructure.setBasePair(left, i, left+mm-1, k);
                            mfeStructure.setStrandRank(t, left + mm - 1);
                            return     backtrackHelper(mm, s, i, t, k, left, left + mm - 1, M(mm, s, i, t, k))
                                    && backtrackHelper(m-mm+1, t, k, r, j+1, left + mm - 1, right, M(m-mm+1, t, k, r, j+1));
                    }}
                }
        }
        System.out.println("Cannot find the follow-up DP decision");
        return false;
    }
    private boolean backtrackHelper(int m, int s, int i, int r, int j, int left, int right, int val){
        if (m == 1){
            if (j - i < 3) return true; // empty region
            if (j > 0 && i+1 < sp.getStrandLength(s) && sp.getM(1, s, i+1, s, j-1) == val){
            return recBacktrack(s, i+1, s, j-1, left, right, val);
        }
        }
        else if (j > 0){
            if (i+1 < sp.getStrandLength(s)){
                if (sp.getM(m, s, i+1, r, j-1) == val){
                    return recBacktrack(s, i+1, r, j-1, left, right, val);
                }
            }
            else if (m == 2){
                if (sp.getM(1, r, 0, r, j-1) == val){
                    return recBacktrack(r, 0, r, j-1, right, right, sp.getM(1, r, 0, r, j-1));
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++) if (sp.getM(m-1, t, 0, r, j-1) == val){
                mfeStructure.setStrandRank(t, left+1);
                return recBacktrack(t, 0, r, j-1, left+1, right, sp.getM(m-1, t, 0, r, j-1));
            }
        }
        else{
            if (i+1 < sp.getStrandLength(s)){
                if (m==2) {
                    if (sp.getM(1, s, i + 1, s, sp.getStrandLength(s) - 1) == val) {
                        return recBacktrack(s, i + 1, s, sp.getStrandLength(s) - 1, left, left, sp.getM(1, s, i + 1, s, sp.getStrandLength(s) - 1));
                    }
                }
                else for (int u = 0; u < sp.getNumStrands(); u++) if (sp.getM(m-1, s, i, u, sp.getStrandLength(u)-1) == val){
                    mfeStructure.setStrandRank(u, right-1);
                    return recBacktrack(s, i, u, sp.getStrandLength(u)-1, left, right-1, sp.getM(m-1, s, i, u, sp.getStrandLength(u)-1));
                }
            }
            if (m <= 2) return true; //empty region
        }
        System.out.printf("!!!!!!! -- %s %s %s %s -- ", left, i, right, j);
        return false;
    }
}
