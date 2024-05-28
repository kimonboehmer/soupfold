# SoupFold

SoupFold allows computing the minimum free energy and the partition function over a set of RNA sequences. Each sequence can appear arbitrarily often in a secondary structure, as long as the total number of strands is $m$. <br>
Furthermore, the code contains the experiments for the paper "RNA Triplet Repeats: Improved Algorithms for1
Structure Prediction and Interactions". <br>

## Using the program

### Creating a strand pool

The input to the algorithm is given as a `StrandPool`. This interface must return the number and length of the strands as well as the letters associated to a position. Furthermore, it stores the DP table and is responsible for initializing, reading and updating it.

#### General pool

The general strand pool is not yet implemented.

#### Triplet pool

You should use this strand pool if all RNA strands are triplet repeats, i.e. of the form $(XYZ)^k$ for $X,Y,Z \in \{A,C,G,U\}$ and $k>0$. The stored DP table uses the repetitive structure of triplet repeats to improve the running time.<br>
There are three constructors:
- `public TripletPool(List<String> strands)`: The input format is a list of strings of tge form $XYZn$, where $X,Y,Z \in \{A,C,G,U\}$ and $n>0$ is a positive integer which denotes the number of repeats of pattern $XYZ$.
- `public TripletPool(Base[] pattern, int mid, int rad)`: This creates a _homogeneous_ triplet pool of pattern `pattern`. Exactly all strand lengths from `mid`-`rad` to `mid`+`rad` are represented in the pool.
- `TripletPool(Base[] pattern, int mid, int rad, int num)`: This is again a homogeneous triplet pool. Compared to the previous constructor, not all strand lengths in the interval from `mid`-`rad` to `mid`+`rad` are represented, but only `num` random ones.

### Creating the DP
To run the algorithm, we just create an instance of `DP`: <br>
`public DP(StrandPool sp, int m, int theta, boolean conn, DPType dpt)` <br>
`m` is the number of interacting strands, `theta`is the minimum base pair span, `conn` should be `true` if we want to enforce the considered structures to be connected and `false` otherwise, and `DPType` specifies whether we want to compute the MFE or the Partition function, and is detailed in the following.<br>
After the constructor is called, we can call two methods:
- `getMFE()` gives the value of the MFE/PF.
- `backtrack` will return a MFE secondary structure in the case that DPType is of type `MFE`, and will sample a secondary structure from the Boltzmann ensemble if `DPType` is of type `PartitionFunction`.

#### MFE

We do not need to give any additional information for the MFE and can just create an instance by <br>
`DPType dpt = new MFE();`

#### Partition Function

For the partition function, we need to specify the temperature in Kelvin. For example: <br>
`DPType dpt = new PartitionFunction(300);`

### Simple example

`LinkedList<String> strands = new LinkedList<>();
strands.add("CAG9");
strands.add("GUU9");
strands.add("ACG9");
StrandPool sp = new TripletPool(strands);
DP mfe = new DP(sp, 3, 3, true, new MFE());
SecondaryStructure st = mfe.backtrack();
System.out.printf("Minimum free energy: %s", mfe.getMFE());
System.out.printf("Secondary structure with minimum free energy: %s", st);
DP pf = new DP(sp, 3, 3, true, new PartitionFunction(300));
SecondaryStructure st = mfe.backtrack();
System.out.printf("Partition Function value: %s", mfe.getMFE());
System.out.printf("Sampled secondary structure from the Boltzmann ensemble: %s", st);`



## Experiments reconstruction

All experiments mentioned in the paper can simply be reconstructed by running the main function. The method <br>
`public static void doExperiment(int num)`<br>
takes one argument `num` which specifies the index between 1 and 7 of the experiment to be reconstructed. Index 0 will execute all experiments in order and is the default option.