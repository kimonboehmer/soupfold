package algorithms;
import datastructures.Base;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;

public class DP {
    private final StrandPool sp;
    public static final int NOT_SET = -100000000;
    private double mfeValue;
    private final DPType dpt;
    private final int startM;
    private final int theta;
    private final boolean conn;
    private SecondaryStructure mfeStructure;

    /**
     * @param sp strand pool on which the DP operates
     * @param m total number of required interacting strands
     * @param theta minimum base pair span
     * @param conn true if connectivity of the secondary structure is required
     * @param dpt MFE or PartitionFunction, depending on the wished computation.
     * Initializes the DP *and* computes the MFE/PF immediately.
     */
    public DP(StrandPool sp, int m, int theta, boolean conn, DPType dpt){
        this.sp = sp;
        this.theta = theta;
        mfeValue = NOT_SET;
        startM = m;
        this.conn = conn;
        this.dpt = dpt;
        mfeValue = computeMFE();
    }

    /**
     * @param m number of remaining strands
     * @param r rightmost strand
     * @param j position on r
     * @return MFE/PF for the region from [t_1,r_j], minimized over all choices for t
     */
    private double minOverStrands(int m, int r, int j){
        double res = dpt.forbidden();
        for (int t = 0; t < sp.getNumStrands();t++){
            double val = dpt.sum(getOrComputeM(m, t, 0, r, j, conn), dpt.strandPenalty(t,sp.getStrandLength(t)));
            res = dpt.min(res, val);
        }
        return res;
    }

    /**
     * @param m number of remaining strand
     * @param s leftmost strand
     * @param i position on s
     * @return MFE/PF for the region from [m,s_i,u_|u|], minimized over all choices for u
     */
    private double minOverSecStrands(int m, int s, int i){
        double res = dpt.forbidden();
        for (int u = 0; u < sp.getNumStrands(); u++){
            double val = dpt.sum(getOrComputeM(m, s, i, u, sp.getStrandLength(u) - 1, conn), dpt.strandPenalty(u,sp.getStrandLength(u)));
            res = dpt.min(res, val);
        }
        return res;
    }

    /**
     * @param m number of remaining strands
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if structure is connected
     * @return MFE/PF for the interval ]m,s_i,r_j,c[.
     */
    private double M(int m, int s, int i, int r, int j, boolean c){
        if (m==0) return dpt.initValue();
        if (m==1) { // if the open region is single-stranded and not empty, compute the included closed region
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
        else{ // structure does not need to be connected and we are at the endpoints of two strands
            double res = dpt.forbidden();
            if (m == 3){ // if we had 3 strands, one central strand is remaining
                for (int t = 0; t < sp.getNumStrands();t++){
                    double val = dpt.sum(getOrComputeM(1, t, 0, t, sp.getStrandLength(t)-1, false), dpt.strandPenalty(t,sp.getStrandLength(t)));
                    res = dpt.min(res, val);
                }
            }
            else{ // else, we have to minimize/sum over two inner strands
                for (int t = 0; t < sp.getNumStrands(); t++) for (int u = 0; u < sp.getNumStrands(); u++){
                    double val = dpt.sum(dpt.sum(getOrComputeM(m-2, t, 0, u, sp.getStrandLength(u) - 1, false), dpt.strandPenalty(t,sp.getStrandLength(t))), dpt.strandPenalty(u,sp.getStrandLength(u)));
                    res = dpt.min(res, val);
                }
            }
            return res;
        }
    }

    /**
     * @param s strand of left base pair endpoint
     * @param i position of left base pair endpoint
     * @param r strand of right base pair endpoint
     * @param j position of right base pair endpoint
     * @return true iff (s_i,r_j) is a valid base pair
     */
    private boolean bp(int s, int i, int r, int j){
        return (Base.pair(sp.getBase(s, i), sp.getBase(r, j)));
    }

    /**
     * @param m number of remaining strands
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if structure is connected
     * @return MFE/PF over the interval [m,s_i,r_j,c] where s_i is paired
     */
    private double computeMultiloop(int m, int s, int i, int r, int j, boolean c){
        double mfe = dpt.forbidden();
        for (int mm = 1; mm <= m; mm++){
            if (mm == 1){ //case 1: interior base pair
                int lim = sp.getStrandLength(s) - 1;
                if (m == 1) lim = j;
                for (int k = i + theta + 1; k <= lim; k++) {
                    if (bp(s,i,s,k)){
                        double val = dpt.sum(dpt.sum(dpt.E(), M(1, s, i, s, k, false)), M(m, s, k, r, j+1, (m > 1) && c));
                        mfe = dpt.min(mfe, val);
                    }
                }
            }
            else if (mm == m){ //case 2: the base pair connects s and r
                for (int k = 0; k <= j; k++) if (bp(s,i,r,k)){
                    double val = dpt.sum(dpt.sum(dpt.E(), M(m, s, i, r, k, false)),M(1, r, k, r, j+1, false));
                    mfe = dpt.min(mfe, val);
                }
            }
            else for (int t = 0; t < sp.getNumStrands(); t++){ //case 3: the base pair connects s and some t
                for (int k = 0; k < sp.getStrandLength(t); k++) if(bp(s, i, t, k)){
                    double val = dpt.sum(dpt.sum(dpt.sum(dpt.E(), M(mm, s, i, t, k, false)), M(m-mm+1, t, k, r, j+1, c)), dpt.strandPenalty(t,sp.getStrandLength(t)));
                    mfe = dpt.min(mfe, val);
                }
            }
        }
        return mfe;
    }

    /**
     * @param m number of remaining strands
     * @param s index of leftmost strand
     * @param i leftmost position on s
     * @param r index of rightmost strand
     * @param j rightmost position on r
     * @param c true if s and r have to be connected, false else
     * @return MFE/PF for interval [m,s_i,r_j,c]
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
        double mfe = dpt.min(unpaired, multi); //minimize/sum over unpaired and paired s_i
        sp.setM(m, s, i, r, j, c, mfe);
        return mfe;
    }

    /**
     * @param m number of remaining strands
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if structure is connected
     * Queries a table entry, and computes it if the entry is not filled.
     * @return MFE/PF for interval [m,s_i,r_j,c]
     */
    private double getOrComputeM(int m, int s, int i, int r, int j, boolean c){
        double val = sp.getM(m, s, i, r, j, c); //query the DP table
        if (val == NOT_SET){ // if the value is not set yet, compute it
            val = computeMFERegion(m, s, i, r, j, c);
        }
        return val;
    }

    /**
     * @return the MFE/PF over the complete set of strands, with the given startM.
     * Fills the DP tables and thus allows backtracking to be called.
     */
    private double computeMFE(){
        sp.initializeTable(startM, theta, dpt.initValue());
        double mfe = dpt.forbidden();
        for (int s = 0; s < sp.getNumStrands(); s++) { // minimize/sum over outermost strands
            if (startM == 1){ // if there is only a single strand, s and r must be equal
                double val = dpt.sum(getOrComputeM(startM, s, 0, s, sp.getStrandLength(s) - 1, false), dpt.strandPenalty(s,sp.getStrandLength(s)));
                mfe = dpt.min(mfe, val);
            }
            else for (int r = 0; r < sp.getNumStrands(); r++) {
                double val = dpt.sum(dpt.sum(getOrComputeM(startM, s, 0, r, sp.getStrandLength(r) - 1, conn), dpt.strandPenalty(s,sp.getStrandLength(s))), dpt.strandPenalty(r,sp.getStrandLength(r)));
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
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(s,sp.getStrandLength(s))))) {
                    mfeStructure.setStrandRank(s, 0);
                    recBacktrack(s, 0, s, sp.getStrandLength(s) - 1,  false, 0, 0, v);
                    return mfeStructure;
                }
            }
        else for (int r = 0; r < sp.getNumStrands(); r++) {
            double v = sp.getM(startM, s, 0, r, sp.getStrandLength(r) - 1, conn);
            if (dpt.btChoose(dpt.sum(dpt.sum(v, dpt.strandPenalty(s,sp.getStrandLength(s))), dpt.strandPenalty(r,sp.getStrandLength(r))))) {
                mfeStructure.setStrandRank(s, 0);
                mfeStructure.setStrandRank(r, startM - 1);
                recBacktrack(s, 0, r, sp.getStrandLength(r) - 1,  conn, 0, startM - 1, v);
                return mfeStructure;
            }
        }
        }
        return mfeStructure;
    }

    /**
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if the structure is connected
     * @param left position in the strand permutation of the leftmost strand
     * @param right position in the strand permutation of the rightmost strand
     * @param val target MFE/PF value.
     * Recursively builds the secondary structure stored in "mfeStructure".
     * @return true if no error occurred
     */
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
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(t,sp.getStrandLength(t))))) {
                    mfeStructure.setStrandRank(t, left + 1);
                    return recBacktrack(t, 0, r, j, conn, left + 1, right, v);
                }
            }
        }
        return backtrackMultiloop(m, s, i, r, j, c, left, right);
    }

    /**
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if the structure is connected
     * @param left position in the strand permutation of the leftmost strand
     * @param right position in the strand permutation of the rightmost strand
     * Recursively builds the secondary structure stored in "mfeStructure",
     * for the case where s_i is paired.
     * val is propagated from the previous call.
     * @return true if no error occurred
     */
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
                        if (bp(s,i,t,k) && dpt.btChoose(dpt.sum(dpt.sum(dpt.sum(dpt.E(),m1), m2), dpt.strandPenalty(t,sp.getStrandLength(t))))){
                            mfeStructure.setBasePair(left, i, left+mm-1, k);
                            mfeStructure.setStrandRank(t, left + mm - 1);
                            return     backtrackHelper(mm, s, i, t, k, false, left, left + mm - 1, m1)
                                    && backtrackHelper(m-mm+1, t, k, r, j+1, c, left + mm - 1, right, m2);
                    }}
                }
        }
        throw new RuntimeException("Inconsistent Backtracking - needs bug fixing!");
    }

    /**
     * @param s leftmost strand
     * @param i position on s
     * @param r rightmost strand
     * @param j position on r
     * @param c true if the structure is connected
     * @param left position in the strand permutation of the leftmost strand
     * @param right position in the strand permutation of the rightmost strand
     * @param val target MFE/PF value.
     * Recursively builds the secondary structure stored in "mfeStructure",
     * for the open interval ]m,s_i,r_j,c[.
     * @return true if no error occurred
     */
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
                if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(t,sp.getStrandLength(t))))) {
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
                    if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(u,sp.getStrandLength(u))))){
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
                    if (dpt.btChoose(dpt.sum(v, dpt.strandPenalty(t,sp.getStrandLength(t))))){
                        mfeStructure.setStrandRank(t, right-1);
                        return recBacktrack(t, 0, t, sp.getStrandLength(t)-1, conn, left+1, right-1, v);
                    }
                }
            }
            for (int t = 0; t < sp.getNumStrands(); t++) for (int u = 0; u < sp.getNumStrands(); u++){
                double v = sp.getM(m-2, t, 0, u, sp.getStrandLength(u)-1, conn);
                if (dpt.btChoose(dpt.sum(dpt.sum(v, dpt.strandPenalty(t,sp.getStrandLength(t))), dpt.strandPenalty(u,sp.getStrandLength(u))))){
                    mfeStructure.setStrandRank(t, left+1);
                    mfeStructure.setStrandRank(u, right-1);
                    return recBacktrack(t, 0, u, sp.getStrandLength(u)-1, conn, left+1, right-1, v);
                }
            }
        }
        throw new RuntimeException("Inconsistent Backtracking - needs bug fixing!!");
    }

    /**
     * @return the MFE or PF value.
     */
    public double getMFE(){
        return mfeValue;
    }
}
