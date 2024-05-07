public interface StrandPool {
int getNumStrands();
int getM(int m, int s, int i, int r, int j);
void setM(int m, int s, int i, int r, int j, int value);
Base getBase(int s, int pos);
int getStrandLength(int s);
void initializeTable(int m, int theta);
String toString(int strand);
}
