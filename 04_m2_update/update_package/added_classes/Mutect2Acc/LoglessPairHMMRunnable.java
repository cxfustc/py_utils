package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.util.concurrent.LinkedBlockingDeque;

public class LoglessPairHMMRunnable implements Runnable {
    private final MonitorRunnable monitor;
    private final boolean isOMPSupported;
    private final M2ArgumentCollection MTAC;
    private final SampleList samplesList;

    private LinkedBlockingDeque<AssemblyResultSet> inputQueue;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;

    LoglessPairHMMRunnable (
            final MonitorRunnable monitor,
            final boolean isOMPSupported,
            final M2ArgumentCollection MTAC,
            final SampleList samplesList
    ) {
        this.monitor = monitor;
        this.isOMPSupported = isOMPSupported;
        this.MTAC = MTAC;
        this.samplesList = samplesList;

        if (isOMPSupported) {
            inputQueue = monitor.loglessPairhMMInputQueue;
        } else {
            inputQueue = monitor.localAssemblyResultQueue;
        }

        createLikelihoodCalculationEngine();
    }

    public void run () {
        monitor.numInitializedThreadsForLoglessPairHMM.incrementAndGet();

        monitor.localAssemblyRunnableLock.lock();
        monitor.numRemainedThreadsForLoglessPairHMM.incrementAndGet();
        monitor.localAssemblyRunnableLock.unlock();

        monitor.runnableIsStarted = true;

        while (true) {
            AssemblyResultSet assemblyResult = inputQueue.poll();
            if (assemblyResult == null) {
                if (monitor.pairHMMSelectionDone) {
                    break;
                }
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
                continue;
            }

            doWork(assemblyResult);

            monitor.numActiveRegionAddressedByPairHMM.incrementAndGet();

            while (!monitor.pairHMMResultQueue.offer(assemblyResult)) {
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
            }

            if (monitor.tooManyEnginesForLoglessPairHMM && monitor.loglessPairHMMRunnableLock.tryLock()) {
                if (monitor.tooManyEnginesForLoglessPairHMM) {
                    monitor.tooManyEnginesForLoglessPairHMM = false;
                    monitor.loglessPairHMMRunnableLock.unlock();
                    break;
                }
                monitor.loglessPairHMMRunnableLock.unlock();
            }
        }

        monitor.loglessPairHMMRunnableLock.lock();
        monitor.numRemainedThreadsForLoglessPairHMM.decrementAndGet();
        monitor.loglessPairHMMRunnableLock.unlock();
    }

    private ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine() {
        final LikelihoodEngineArgumentCollection likelihoodArgs = MTAC.likelihoodArgs;
        final double log10GlobalReadMismappingRate = likelihoodArgs.phredScaledGlobalReadMismappingRate<0 ? -Double.MAX_VALUE
                : QualityUtils.qualToErrorProbLog10(likelihoodArgs.phredScaledGlobalReadMismappingRate);

        switch (likelihoodArgs.likelihoodEngineImplementation) {
            case PairHMM:
                likelihoodCalculationEngine = new PairHMMLikelihoodCalculationEngine(
                        (byte)likelihoodArgs.gcpHMM,
                        likelihoodArgs.pairHMMNativeArgs.getPairHMMArgs(),
                        likelihoodArgs.pairHMM,
                        log10GlobalReadMismappingRate,
                        likelihoodArgs.pcrErrorModel,
                        likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD,
                        isOMPSupported);
                break;
            case Random:
                likelihoodCalculationEngine = new RandomLikelihoodCalculationEngine();
                break;
            default:
                throw new UserException("Unsupported likelihood calculation engine.");
        }

        return likelihoodCalculationEngine;
    }

    private void doWork (final AssemblyResultSet assemblyResult) {
        assemblyResult.readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, assemblyResult.reads);
    }
}
