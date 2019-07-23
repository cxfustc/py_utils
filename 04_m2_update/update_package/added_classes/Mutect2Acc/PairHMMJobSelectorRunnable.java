package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;

public class PairHMMJobSelectorRunnable implements Runnable {
    private final MonitorRunnable monitor;

    private final int COMPUTATIONAL_INTENSIVE_LIMIT = 800;
    private final int MIN_COMPUTATIONAL_INTENSIVE_LIMIT = 200;
    private final int ADJUST_STEP = 200;

    PairHMMJobSelectorRunnable (
            MonitorRunnable monitor
    ) {
        this.monitor = monitor;
    }

    public void run () {
        while (true) {
            AssemblyResultSet assemblyResult = monitor.localAssemblyResultQueue.poll();
            if (assemblyResult == null) {
                if (monitor.numRemainedThreadsForLocalAssembly.get() <= 0) {
                    break;
                }
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
                continue;
            }

            if (isComputationallyIntensive(assemblyResult)) {
                monitor.numActiveRegionAssignedToVector.incrementAndGet();
                while (!monitor.vectorLoglessPairHMMInputQueue.offer(assemblyResult)) {
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {

                    }
                }
            } else {
                monitor.numActiveRegionAssignedToLogless.incrementAndGet();
                while (!monitor.loglessPairhMMInputQueue.offer(assemblyResult)) {
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {

                    }
                }
            }
        }
        monitor.pairHMMSelectionDone = true;
    }

    private boolean isComputationallyIntensive (final AssemblyResultSet assemblyResult) {
        int readCount = assemblyResult.getRegionForGenotyping().getReads().size();
        int haplotypeCount = assemblyResult.getHaplotypeCount();
        if (readCount*haplotypeCount > COMPUTATIONAL_INTENSIVE_LIMIT) {
            return true;
        } else {
            return false;
        }
    }
}
