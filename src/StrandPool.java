public interface StrandPool {
int getNumStrands();
double getM(int m, int s, int i, int r, int j, boolean c, int diff);
void setM(int m, int s, int i, int r, int j, boolean c, int diff, double value);
Base getBase(int s, int pos);
int getStrandLength(int s);
void initializeTable(int m, int theta, double initValue, int avg);
String toString(int strand);
int getMinDiff();
int getMaxDiff();

}
