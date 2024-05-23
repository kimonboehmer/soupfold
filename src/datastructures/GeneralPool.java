package datastructures;

import datastructures.Base;

public class GeneralPool {
    Base[][] pool;
    int[] appearances;
    public GeneralPool(Base[][] pool, int[] appearances){
        this.pool = pool;
        this.appearances = appearances;
    }
    public GeneralPool(Base[] bases, int mid, int rad, int num){
        assert mid > rad;
        pool = new Base[2 * rad + 1][];
        appearances = new int[2 * rad + 1];
        pool[0] = fillStrandWithTriplet(bases, mid);
        appearances[0] = num;
        for (int i = 1; i <= rad; i++){
            pool[2 * i - 1] = fillStrandWithTriplet(bases, mid - i);
            pool[2 * i] = fillStrandWithTriplet(bases, mid + i);
            appearances[2 * i - 1] = num;
            appearances[2 * i] = num;
        }
    }
    private Base[] fillStrandWithTriplet(Base[] bases, int repeats){
        Base[] triplet = new Base[repeats * bases.length];
        for (int i = 0; i < repeats; i++){
            System.arraycopy(bases, 0, triplet, bases.length * i, bases.length);
        }
        return triplet;
    }
}
