package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

class Mutect2AccEngine {
    public ExecutorService threadPool;
    public Runnable runnable;
    public Future<?> future;

    public Mutect2AccEngine (
            final ExecutorService threadPool,
            final Runnable runnable,
            final Future<?> future
    ) {
        this.threadPool = threadPool;
        this.runnable = runnable;
        this.future = future;
    }
}
