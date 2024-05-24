package experiments;

import algorithms.DP;
import algorithms.PartitionFunction;
import datastructures.SecondaryStructure;
import datastructures.StrandPool;
import datastructures.TripletPool;

import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;

public class Sampling {
    public static SecondaryStructure sampleAndReject(DP dp){
        Random r = new Random();
        SecondaryStructure s;
        LinkedList<Integer>[] connComps;
        double threshold;
        do{
            s = dp.backtrack();
            connComps = connectedComponents(s);
            threshold = 1.0 / (double) computeOvercountingFactor(connComps);
        }
        while (/*connComps.length > 1 && */r.nextDouble() > threshold);
        return s;
    }
    private static LinkedList<Integer>[] connectedComponents(SecondaryStructure st){
        LinkedList<LinkedList<Integer>> connComps = new LinkedList<>();
        int m = st.getM();
        boolean[] visited = new boolean[m];
        for (int s = 0; s < m; s++) if (!visited[s]){
            LinkedList<Integer> connComp = new LinkedList<>();
            Queue<Integer> q = new LinkedList<>();
            q.offer(s);
            visited[s] = true;
            while (!q.isEmpty()){
                int r = q.poll();
                connComp.add(r);
                for (int i = 0; i < st.getStrandPool().getStrandLength(st.getStrandFromPosition(r)); i++){
                    int p = st.getPairedStrand(r, i);
                    if (p >= 0 && !visited[p]){
                        q.offer(p);
                        visited[p] = true;
                    }
                }
            }
            connComps.add(connComp);
        }
        LinkedList<Integer>[] arr = new LinkedList[connComps.size()];
        int i = 0;
        for (LinkedList<Integer> connComp : connComps) arr[i++] = connComp;
        return arr;
    }
    private static int computeOvercountingFactor(LinkedList<Integer>[] connComps){
        int prod = 1;
        for (int i = 0; i < connComps.length; i++){
            int sum = 0;
            for (int j = 0; j < i; j++) sum += connComps[j].size();
            if (sum == 0) sum = 1;
            prod *= sum * connComps[i].size();
        }
        return prod;
    }
    public static double[][][][] basePairProbabilities(StrandPool sp, int m, int sampleSize){
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        int maxStrandLength = 0;
        for (int i = 0; i < sp.getNumStrands(); i++) maxStrandLength = Math.max(maxStrandLength, sp.getStrandLength(i));
        double[][][][] table = new double[sp.getNumStrands()][maxStrandLength][sp.getNumStrands()][maxStrandLength];
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure st = dp.backtrack();
            for (int mm = 0; mm < m; mm++){
                int s = st.getStrandFromPosition(mm);
                for (int i = 0; i < sp.getStrandLength(s); i++){
                    int r = st.getPairedStrand(mm, i);
                    if (r >= 0){
                        table[s][i][st.getStrandFromPosition(r)][st.getPairedIndex(mm, i)]++;
                    }
                }
            }
        }
        for (int s = 0; s < sp.getNumStrands(); s++) for (int i = 0; i < sp.getStrandLength(s); i++) for (int r = 0; r < sp.getNumStrands(); r++) for (int j = 0; j < sp.getStrandLength(r); j++){
            table[s][i][r][j] /= sampleSize;
        }
        return table;
    }
    public static double [][] bpTypes(StrandPool sp, int m, int sampleSize){
        double[][] data = new double[sp.getNumStrands()][sp.getNumStrands()];
        double[][][][] table = basePairProbabilities(sp, m, sampleSize);
        for (int s = 0; s < sp.getNumStrands();s++) for (int r = 0; r < sp.getNumStrands(); r++){
            double sum = 0;
            for (int i = 0; i < sp.getStrandLength(s);i++) for (int j = 0; j < sp.getStrandLength(r);j++){
                sum += table[s][i][r][j];
            }
            if (s == r) sum /= 2;
            data[s][r] = sum;
            System.out.printf("From strand %s to strand %s: %f base pairs in average.\n", s, r, sum);
        }
        return data;
    }
    public static double[][][][] basePairProbabilitiesByPosition(StrandPool sp, int m, int sampleSize){
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        int maxStrandLength = 0;
        for (int i = 0; i < sp.getNumStrands(); i++) maxStrandLength = Math.max(maxStrandLength, sp.getStrandLength(i));
        double[][][][] table = new double[m][maxStrandLength][m][maxStrandLength];
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure st = dp.backtrack();
            for (int mm = 0; mm < m; mm++) for (int i = 0; i < sp.getStrandLength(st.getStrandFromPosition(mm)); i++){
                int r = st.getPairedStrand(mm, i);
                if (r >= 0){
                    table[mm][i][r][st.getPairedIndex(mm, i)]++;
                }
            }
        }
        for (int s = 0; s < m; s++) for (int i = 0; i < table[s].length; i++) for (int r = 0; r < m; r++) for (int j = 0; j < table[s][i][r].length; j++){
            table[s][i][r][j] /= sampleSize;
        }
        return table;
    }
    public static double interactionProbability(StrandPool sp, int m, int sampleSize){
        double sum = 0;
        double homoSum = 0;
        int norm = 0;
        double[][][][] table = basePairProbabilitiesByPosition(sp, m, sampleSize);
        for (int s = 0; s < m; s++) for (int i = 0; i < table[s].length; i++){
            double countAll = 0;
            double countInteriorBPs = 0;
            double countHomoBPs = 0;
            for (int r = 0; r < m; r++) for (int j = 0; j < table[s][i][r].length; j++){
                countAll += table[s][i][r][j];
            }
            if (countAll > 0){
                for (int j = 0; j < table[s][i][s].length; j++) countInteriorBPs += table[s][i][s][j];
                sum += countInteriorBPs / countAll;
                norm++;
            }
        }
        return 1 - (sum / (double) norm);
    }
    public static double[] expNumOccurencesOfStrands(StrandPool sp, int m, int sampleSize){
        double[] data = new double[sp.getNumStrands()];
        DP dp = new DP(sp, m, 3, false, new PartitionFunction(300));
        for (int l = 0; l < sampleSize; l++) {
            SecondaryStructure st = dp.backtrack();
            for (int i = 0; i < m; i++) data[st.getStrandFromPosition(i)]++;
        }
        for (int i = 0; i < sp.getNumStrands(); i++){
            data[i] /= sampleSize;
            System.out.println(data[i]);
        }
        return data;
    }
    public static LinkedList<Double> connectivityExperiment(StrandPool sp, int m, int sampleSize){
        int[] ccSize = new int[m+1];
        DP dp  = new DP(sp, m, 3, false, new PartitionFunction(1000));
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure s = sampleAndReject(dp);
            LinkedList<Integer>[] cc = connectedComponents(s);
            ccSize[cc.length]++;
        }
        LinkedList<Double> ll = new LinkedList<>();
        for (int i = 0; i < ccSize.length; i++){
            System.out.printf("num of cc = %d: %d times\n", i, ccSize[i]);
            ll.add((double) ccSize[i]);
        }
        ll.set(0, 1 - ll.get(1) / (double) sampleSize);
        return ll;
    }
    public static double avgBPs(StrandPool sp, int m, int sampleSize){
        double sum = 0;
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure s = dp.backtrack();
            if (l == 5000) System.out.println(s.toString());
            int len = 0;
            for (int i = 0; i < m; i++) len += sp.getStrandLength(s.getStrandFromPosition(i));
            sum += s.getNumBPs() / (double) len;
        }
        return sum / (double) sampleSize;
    }
    public static double[] classifyBasePairs(TripletPool sp, int m, int sampleSize){
        int interior = 0;
        int homoExt = 0;
        int heteroExt = 0;
        DP dp = new DP(sp, m, 3, true, new PartitionFunction(300));
        for (int l = 0; l < sampleSize; l++){
            SecondaryStructure st = dp.backtrack();
            for (int s = 0; s < m; s++){
                for (int i = 0; i < sp.getStrandLength(st.getStrandFromPosition(s)); i++){
                    int r = st.getPairedStrand(s, i);
                    if (r < 0) continue;
                    if (s == r) interior++;
                    else if (sp.getPattern(st.getStrandFromPosition(s)) == sp.getPattern(st.getStrandFromPosition(r))) homoExt++;
                    else heteroExt++;
                }
            }
        }
        double sum = interior + homoExt + heteroExt;
        return new double[]{interior / sum, homoExt / sum, heteroExt / sum};
    }
}
