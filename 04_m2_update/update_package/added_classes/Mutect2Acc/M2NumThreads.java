package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

public class M2NumThreads {
    public int numThreadsForDataPrepare;
    public int numThreadsForLocalAssembly;
    public int numThreadsForLoglessPairHMM;
    public int numThreadsForGenotyping;

    public int minNumThreadsForDataPrepare;
    public int minNumThreadsForLocalAssembly;
    public int minNumThreadsForLoglessPairHMM;
    public int minNumThreadsForGenotyping;

    public int maxNumThreadsForDataPrepare;
    public int maxNumThreadsForLocalAssembly;
    public int maxNumThreadsForLoglessPairHMM;
    public int maxNumThreadsForGenotyping;
}
