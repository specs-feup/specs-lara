package pt.up.fe.specs.lara.kadabra.utils;


public class CountingMonitorList {
    private int[] counter;

    public CountingMonitorList(int size) {
        counter = new int[size];
    }

    public void increment(int id) {
        counter[id]++;
    }

    public int get(int id) {
        return counter[id];
    }

    public void reset(int id) {
        counter[id] = 0;
    }
}

