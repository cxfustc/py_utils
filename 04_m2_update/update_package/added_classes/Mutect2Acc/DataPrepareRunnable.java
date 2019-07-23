package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;

import java.util.Iterator;
import java.util.List;

public class DataPrepareRunnable implements Runnable {
    private ReferenceDataSource reference;
    private ReadsDataSource reads;
    private FeatureManager features;

    private ReadsDownsampler readsDownsampler;
    private SAMFileHeader header;

    private List<MultiIntervalLocalReadShard> readShards;
    private int numShards;

    private MonitorRunnable monitor;
    private CountingReadFilter countedFilter;
    private M2IsActiveEngine engine;

    private final ReferenceInputArgumentCollection referenceArguments;
    private final M2ArgumentCollection MTAC;
    private final int minRegionSize;
    private final int maxRegionSize;
    private final int assemblyRegionPadding;
    private final int maxProbPropagationDistance;
    private final double activeProbThreshold;
    private final boolean includeReadsWithDeletionsInIsActivePileups;

    DataPrepareRunnable (
            final List<MultiIntervalLocalReadShard> readShards,
            final CommandLineProgram instance,
            final MonitorRunnable monitor,
            final CountingReadFilter countedFilter,
            final ReferenceInputArgumentCollection referenceArguments,
            final ReadInputArgumentCollection readArguments,
            final M2ArgumentCollection MTAC,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer,
            final int maxReadsPerAlignmentStart,
            final Logger logger,
            final String tumorSample,
            final String normalSample,
            final int minRegionSize,
            final int maxRegionSize,
            final int assemblyRegionPadding,
            final double activeProbThreshold,
            final int maxProbPropagationDistance,
            final boolean includeReadsWithDeletionsInIsActivePileups
    ) {
       reference = Mutect2Acc.LocalInitializeReference(referenceArguments);
       reads = Mutect2Acc.LocalInitializeReads(referenceArguments,readArguments,cloudPrefetchBuffer,cloudIndexPrefetchBuffer);
       features = Mutect2Acc.LocalInitializeFeatures(instance, referenceArguments, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);

       this.readShards = readShards;
       this.monitor = monitor;
       this.countedFilter = countedFilter;

       this.referenceArguments = referenceArguments;
       this.MTAC = MTAC;
       this.minRegionSize = minRegionSize;
       this.maxRegionSize = maxRegionSize;
       this.assemblyRegionPadding = assemblyRegionPadding;
       this.maxProbPropagationDistance = maxProbPropagationDistance;
       this.activeProbThreshold = activeProbThreshold;
       this.includeReadsWithDeletionsInIsActivePileups = includeReadsWithDeletionsInIsActivePileups;

       numShards = readShards.size();
       readsDownsampler = new MutectDownsampler(maxReadsPerAlignmentStart, MTAC.maxSuspiciousReadsPerAlignmentStart, MTAC.downsamplingStride);

       header = reads.getHeader();
       engine = new M2IsActiveEngine(MTAC, header, logger, tumorSample, normalSample);
    }

    public void run () {
        monitor.numInitializedThreadsForDataPrepare.incrementAndGet();

        monitor.dataPrepareRunnableLock.lock();
        monitor.numRemainedThreadsForDataPrepare.incrementAndGet();
        monitor.dataPrepareRunnableLock.unlock();

        monitor.runnableIsStarted = true;

        while (true) {
            int currentIndex = monitor.indexOfCurrentShard.getAndIncrement();
            if (currentIndex >= numShards) {
                break;
            }

            MultiIntervalLocalReadShard readShard = readShards.get(currentIndex);

            readShard.setPreReadFilterTransformer(ReadTransformer.identity());
            readShard.setReadsSource(reads);
            readShard.setReadFilter(countedFilter);
            readShard.setDownsampler(readsDownsampler);
            readShard.setPostReadFilterTransformer(ReadTransformer.identity().andThen(Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), !MTAC.dontClipITRArtifacts)));

            Iterator<AssemblyRegion> assemblyRegionIterator = new AssemblyRegionIterator(
                    readShard,
                    header,
                    reference,
                    features,
                    engine,
                    minRegionSize,
                    maxRegionSize,
                    assemblyRegionPadding,
                    activeProbThreshold,
                    maxProbPropagationDistance,
                    includeReadsWithDeletionsInIsActivePileups
            );

            int indexOfAssemblyRegion = 0;
            while (assemblyRegionIterator.hasNext()) {
                AssemblyRegion assemblyRegion = assemblyRegionIterator.next();

                if (!assemblyRegion.isActive() || assemblyRegion.size()==0) {
                    continue;
                }

                //String s = String.format("%s:%d-%d\tread count: %d",
                //        assemblyRegion.getContig(), assemblyRegion.getStart(), assemblyRegion.getEnd(),
                //        assemblyRegion.getReads().size());
                //System.err.println (s);

                monitor.numActiveRegionCreatedByDataPrepare.incrementAndGet();
                while (!monitor.activeRegionQueue.offer(assemblyRegion)) {
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {

                    }
                }
            }
            assemblyRegionIterator = null;
            monitor.numActiveRegionOnShards[currentIndex] = indexOfAssemblyRegion;

            monitor.numRemainedIntervals.decrementAndGet();

            if (monitor.tooManyEnginesForDataPrepare && monitor.dataPrepareRunnableLock.tryLock()) {
                if (monitor.tooManyEnginesForDataPrepare) {
                    monitor.tooManyEnginesForDataPrepare = false;
                    monitor.dataPrepareRunnableLock.unlock();
                    break;
                }
                monitor.dataPrepareRunnableLock.unlock();
            }
        }

        monitor.dataPrepareRunnableLock.lock();
        monitor.numRemainedThreadsForDataPrepare.decrementAndGet();
        monitor.dataPrepareRunnableLock.unlock();
    }
}
