package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.AssemblyRegionWalker;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc.Mutect2Acc;

/**
 * Arguments for native PairHMM implementations
 */
public class PairHMMNativeArgumentCollection {

    @Argument(fullName = "native-pair-hmm-use-double-precision", doc="use double precision in the native pairHmm. " +
            "This is slower but matches the java implementation better", optional = true)
    private boolean useDoublePrecision = false;

    public PairHMMNativeArguments getPairHMMArgs(){
        final PairHMMNativeArguments args = new PairHMMNativeArguments();
        if (Mutect2Acc.pairHmmNativeThreads <= 0) {
            Mutect2Acc.pairHmmNativeThreads = Runtime.getRuntime().availableProcessors() / 2;
        }
        args.maxNumberOfThreads = Mutect2Acc.pairHmmNativeThreads;
        args.useDoublePrecision = useDoublePrecision;
        return args;
    }

}
