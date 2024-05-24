import experiments.Helper;

import java.io.IOException;

import static experiments.FinalExperiments.doExperiment;

public class Main {
    public static void main(String[] args) throws IOException {
        for (int m = 2; m < 6; m++) for (int t = 0; t < 3; t++) {
            double[][] data = Helper.readCSV("table_m" + m + "_t" + t + ".csv");
            data = Helper.changeZeros(data);
            Helper.writeCSV(data, "table_m" + m + "_t" + t + "_new.csv");
        }
    }
}
