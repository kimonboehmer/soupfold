package algorithms;
import datastructures.Base;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;

public class DP {
    private final StrandPool sp;
    public static final int NOT_SET = 100000;
    private double mfeValue;
    private final DPType dpt;
    private final int startM;
    private final int theta;
    private final boolean conn;
    private SecondaryStructure mfeStructure;
    public DP(StrandPool sp, int m, int theta, boolean conn, DPType dpt){
        this.sp = sp;
        this.theta = theta;
        mfeValue = NOT_SET;
        startM = m;
        this.conn = conn;
        this.dpt = dpt;
        mfeValue = computeMFE();
    }
    private double minOverStrands(int m, int r, int j){
        double res = dpt.forbidden();
        for (int t = 0; t < sp.getNumStrands();t++){
            double val = dpt.sum(getOrComputeM(m, t, 0, r, j, conn), dpt.strandPenalty(sp.getStrandLength(t)));
            res = dpt.min(res, val);
        }
        return res;
    }
    private double minOverSecStrands(int m, int s, int i){
        double res = dpt.forbidden();
        for (int u = 0; u < sp.getNumStrands(); u++){
            double val = dpt.sum(getOrComputeM(m, s, i, u, sp.getStrandLength(u) - 1, conn), dpt.strandPenalty(sp.getStrandLength(u)));
            res = dpt.min(res, val);
        }
        return res;
    }
    private double M(int m, int s, int i, int r, int j, boolean c){
        if (m==0) return dpt.initValue();
        if (m==1) {
            if (i+1 < sp.getStrandLength(s) && j > 0 && j > i+1) return getOrComputeM(1, s, i+1, s, j-1, false);
            return dpt.initValue();
        }
        if (j > 0){
            if (i+1 < sp.getStrandLength(s)) return getOrComputeM(m, s, i+1, r, j-1, c);
            // we leave strand s. If connectivity was required, this is not allowed
            if (c) return dpt.forbidden();
            if (m==2) return getOrComputeM(1, r, 0, r, j-1, false);
            return minOverStrands(m-1, r, j-1);
        }
        // we leave strand r
        if (c) return dpt.forbidden();
        if (i+1 < sp.getStrandLength(s)){
            if (m==2) return getOrComputeM(1, s, i+1, s, sp.getStrandLength(s)-1, false);
            return minOverSecStrands(m-1, s, i+1);
        }
        if (m == 2) return dpt.initValue();
        if (conn) return dpt.forbidden();
        else{
            double res = dpt.forbidden();
            if (m == 3){
                for (int t = 0; t < sp.getNumStrands();t++){
                    double val = dpt.sum(getOrComputeM(1, t, 0, t, sp.getStrandLength(t)-1, false), dpt.strandPenalty(sp.getStrandLength(t)));
                    res = dpt.min(res, val);
                }
            }
            else{
                for (int t = 0; t < sp.getNumStrands(); t++) for (int u = 0; u < sp.getNumStrands(); u++){
                    double val = dpt.sum(dpt.sum(getOrComputeM(m-2, t, 0, u, sp.getStrandLength(u) - 1, false), dpt.strandPenalty(sp.getStrandLength(t))), dpt.strandPenalty(sp.getStrandLength(u)));
                    res = dpt.min(res, val);
                }
            }
            return res;
        }
    }
    private boolean bp(int s, int i, int r, int j){
        return (Base.pair(sp.getBase(s, i), sp.getBase(r, j)));
    }
    private double computeMultiloop(int m, int s, int i, int r, int j, boolean c){
        double mfe = dpt.forbidden();
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s) - 1;
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k <= lim; k++) {
                    if (bp(s,i,s,k)){
                        double val = dpt.sum(dpt.sum(dpt.E(), M(1, s, i, s, k, false)), M(m, s, k, r, j+1, (m > 1) && c));
                        mfe = dpt.min(mfe, val);
                    }
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++) if (bp(s,i,r,k)){
                    double val = dpt.sum(dpt.sum(dpt.E(), M(m, s, i, r, k, false)),M(1, r, k, r, j+1, false));
                    mfe = dpt.min(mfe, val);
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                for (int k = 0; k < sp.getStrandLength(t); k++) if(bp(s, i, t, k)){
                    double val = dpt.sum(dpt.sum(dpt.sum(dpt.E(), M(mm, s, i, t, k, false)), M(m-mm+1, t, k, r, j+1, c)), dpt.strandPenalty(sp.getStrandLength(t)));
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
        else {
            if (c) unpaired = dpt.forbidden();
            else if (m==2) unpaired = getOrComputeM(1, r, 0, r, j, false);
            else unpaired = minOverStrands(m - 1, r, j);
        }
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
    private double computeMFE(){
        sp.initializeTable(startM, theta, dpt.initValue());
        double mfe = dpt.forbidden();
        for (int s = 0; s < sp.getNumStrands(); s++) {
            if (startM == 1){
                double val = dpt.sum(getOrComputeM(startM, s, 0, s, sp.getStrandLength(s) - 1, false), dpt.strandPenalty(sp.getStrandLength(s)));
                mfe = dpt.min(mfe, val);
            }
            else for (int r = 0; r < sp.getNumStrands(); r++) {
                double val = dpt.sum(dpt.sum(getOrComputeM(startM, s, 0, r, sp.getStrandLength(r) - 1, conn), dpt.strandPenalty(sp.getStrandLength(s))), dpt.strandPenalty(sp.getStrandLength(r)));
                mfe = dpt.min(mfe, val);
            }
        }
        mfeValue = mfe;
        return mfe;
    }

    /**
     * @return one optimal secondary structure if dpt is of type MFE, and
     *         one sampled secondary structure if dpt is of type PartitionFunction.
     */
    public SecondaryStructure backtrack(){
        mfeStructure = new SecondaryStructure(sp, startM);
        dpt.btInit(mfeValue);
        for (int s = 0; s < sp.getNumStrands(); s++){
            if (startM == 1){ // in case m = 1, the start- and end-strand has to be the same
                double v = sp.getM(startM, s, 0, s, sp.getStrandLength(s) - 1, false);
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(s))))) {
                    mfeStructure.setStrandRank(s, 0);
                    recBacktrack(s, 0, s, sp.getStrandLength(s) - 1,  false, 0, 0, v);
                    return mfeStructure;
                }
            }
        else for (int r = 0; r < sp.getNumStrands(); r++) {
            double v = sp.getM(startM, s, 0, r, sp.getStrandLength(r) - 1, conn);
            if (dpt.btChoose(dpt.sum(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(s))), dpt.strandPenalty(sp.getStrandLength(r))))) {
                mfeStructure.setStrandRank(s, 0);
                mfeStructure.setStrandRank(r, startM - 1);
                recBacktrack(s, 0, r, sp.getStrandLength(r) - 1,  conn, 0, startM - 1, v);
                return mfeStructure;
            }
        }
        }
        return mfeStructure;
    }
    private boolean recBacktrack(int s, int i, int r, int j, boolean c, int left, int right, double val) {
        dpt.btInit(val);
        if (val == 0) return true;
        if (left == right && j - i <= theta) return true;
        int m = right - left + 1;
        if (val != sp.getM(m, s, i, r, j, c)) throw new RuntimeException("Inconsistent Backtracking - needs bug fixing.");
        // unpaired
        if (i + 1 < sp.getStrandLength(s)){
            double v = sp.getM(m, s, i + 1, r, j, c);
            if (dpt.btChoose(v)) return recBacktrack(s, i + 1, r, j, c, left, right, v);
        }
        else if (!c){
            if (m == 2) {
                double v = sp.getM(1, r, 0, r, j, false);
                if (dpt.btChoose(v)) return recBacktrack(r, 0, r, j, false, right, right, v);
            }
            else for (int t = 0; t < sp.getNumStrands(); t++) {
                double v = sp.getM(m - 1, t, 0, r, j, conn) ;
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(t))))) {
                    mfeStructure.setStrandRank(t, left + 1);
                    return recBacktrack(t, 0, r, j, conn, left + 1, right, v);
                }
            }
        }
        return backtrackMultiloop(m, s, i, r, j, c, left, right);
    }
    private boolean backtrackMultiloop(int m, int s, int i, int r, int j, boolean c, int left, int right){
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){
                int lim = sp.getStrandLength(s) - 1;
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k <= lim; k++){
                    double m1 = M(1, s, i, s, k, false);
                    double m2 = M(m, s, k, r, j+1, m>1 && c);
                    if (bp(s,i,s,k) && dpt.btChoose(dpt.sum(dpt.sum(dpt.E(), m1), m2))){
                        mfeStructure.setBasePair(left, i, left, k);
                        return     backtrackHelper(1, s, i, s, k, false, left, left, m1)
                                && backtrackHelper(m, s, k, r, j+1, (m>1 && c), left, right, m2);
                }
                }
            }
            else if (mm == m){
                for (int k = 0; k <= j; k++){
                    double m1 = M(m, s, i, r, k, false);
                    double m2 = M(1, r, k, r, j+1, false);
                    if (bp(s,i,r,k) && dpt.btChoose(dpt.sum(dpt.sum(dpt.E(), m1),m2))){
                        mfeStructure.setBasePair(left, i, right, k);
                        return     backtrackHelper(m, s, i, r, k, false, left, right, m1)
                                && backtrackHelper(1, r, k, r, j+1, false, right, right, m2);
                }}
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){
                    for (int k = 0; k < sp.getStrandLength(t); k++){
                        double m1 = M(mm, s, i, t, k, false);
                        double m2 = M(m-mm+1, t, k, r, j+1, c);
                        if (bp(s,i,t,k) && dpt.btChoose(dpt.sum(dpt.sum(dpt.sum(dpt.E(),m1), m2), dpt.strandPenalty(sp.getStrandLength(t))))){
                            mfeStructure.setBasePair(left, i, left+mm-1, k);
                            mfeStructure.setStrandRank(t, left + mm - 1);
                            return     backtrackHelper(mm, s, i, t, k, false, left, left + mm - 1, m1)
                                    && backtrackHelper(m-mm+1, t, k, r, j+1, c, left + mm - 1, right, m2);
                    }}
                }
        }
        throw new RuntimeException("Inconsistent Backtracking - needs bug fixing!");
    }
    private boolean backtrackHelper(int m, int s, int i, int r, int j, boolean c, int left, int right, double val){
        dpt.btInit(val);
        if (m == 0) return true;
        if (m == 1){
            if (j - i < theta) return true; // empty region
            if (j > 0 && i+1 < sp.getStrandLength(s)){
                double v = sp.getM(1, s, i+1, s, j-1, false);
                if (dpt.btChoose(v)){
                    return recBacktrack(s, i+1, s, j-1, false, left, right, v);
                }
        }
        }
        else if (j > 0) {
            if (i + 1 < sp.getStrandLength(s)) {
                double v = sp.getM(m, s, i + 1, r, j - 1, c);
                if (dpt.btChoose(v)) {
                    return recBacktrack(s, i + 1, r, j - 1, c, left, right, v);
                }
            } else if (!c){
                if (m == 2) {
                double v = sp.getM(1, r, 0, r, j - 1, false);
                if (dpt.btChoose(v)) {
                    return recBacktrack(r, 0, r, j - 1, false, right, right, v);
                }
            } else for (int t = 0; t < sp.getNumStrands(); t++) {
                double v = sp.getM(m - 1, t, 0, r, j - 1, conn);
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(t))))) {
                    mfeStructure.setStrandRank(t, left + 1);
                    return recBacktrack(t, 0, r, j - 1, conn, left + 1, right, v);
                }
            }
        }
        }
        else if (!c){
            if (i+1 < sp.getStrandLength(s)){
                if (m==2) {
                    double v = sp.getM(1, s, i + 1, s, sp.getStrandLength(s) - 1, false);
                    if (dpt.btChoose(v)) {
                        return recBacktrack(s, i + 1, s, sp.getStrandLength(s) - 1, false, left, left, v);
                    }
                }
                else for (int u = 0; u < sp.getNumStrands(); u++){
                    double v = sp.getM(m-1, s, i+1, u, sp.getStrandLength(u)-1, conn);
                    if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(u))))){
                        mfeStructure.setStrandRank(u, right-1);
                        return recBacktrack(s, i+1, u, sp.getStrandLength(u)-1, conn, left, right-1, v);
                    }
                }
            }
            if (m <= 2) return true; //empty region
            if (conn) throw new RuntimeException("hm");
            if (m == 3){
                for (int t = 0; t < sp.getNumStrands(); t++){
                    double v = sp.getM(1, t, 0, t, sp.getStrandLength(t)-1, conn);
                    if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(t))))){
                        mfeStructure.setStrandRank(t, right-1);
                        return recBacktrack(t, 0, t, sp.getStrandLength(t)-1, conn, left+1, right-1, v);
                    }
                }
            }
            for (int t = 0; t < sp.getNumStrands(); t++) for (int u = 0; u < sp.getNumStrands(); u++){
                double v = sp.getM(m-2, t, 0, u, sp.getStrandLength(u)-1, conn);
                if (dpt.btChoose(dpt.sum(dpt.sum(v, dpt.strandPenalty(sp.getStrandLength(t))), dpt.strandPenalty(sp.getStrandLength(u))))){
                    mfeStructure.setStrandRank(t, left+1);
                    mfeStructure.setStrandRank(u, right-1);
                    return recBacktrack(t, 0, u, sp.getStrandLength(u)-1, conn, left+1, right-1, v);
                }
            }
        }
        throw new RuntimeException("Inconsistent Backtracking - needs bug fixing!!");
    }
    public double getMFE(){
        return mfeValue;
    }
}
