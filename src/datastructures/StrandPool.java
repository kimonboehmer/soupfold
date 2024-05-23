package datastructures;
public interface StrandPool {
int getNumStrands();
double getM(int m, int s, int i, int r, int j, boolean c);
void setM(int m, int s, int i, int r, int j, boolean c, double value);
Base getBase(int s, int pos);
int getStrandLength(int s);
void initializeTable(int m, int theta, double initValue);
String toString(int strand);

}
