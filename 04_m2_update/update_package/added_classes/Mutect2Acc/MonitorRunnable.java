package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.MultiIntervalLocalReadShard;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import scala.tools.nsc.backend.icode.Primitives;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class MonitorRunnable implements Runnable {
    protected AtomicInteger numRemainedThreadsForDataPrepare;
    protected AtomicInteger numRemainedThreadsForLocalAssembly;
    protected AtomicInteger numRemainedThreadsForLoglessPairHMM;
    protected AtomicInteger numRemainedThreadsForGenotyping;

    protected AtomicInteger numInitializedThreadsForDataPrepare;
    protected AtomicInteger numInitializedThreadsForLocalAssembly;
    protected AtomicInteger numInitializedThreadsForLoglessPairHMM;
    protected AtomicInteger numInitializedThreadsForGenotyping;

    protected LinkedBlockingDeque<AssemblyRegion> activeRegionQueue;
    protected LinkedBlockingDeque<AssemblyResultSet> localAssemblyResultQueue;
    protected LinkedBlockingDeque<AssemblyResultSet> loglessPairhMMInputQueue;
    protected LinkedBlockingDeque<AssemblyResultSet> vectorLoglessPairHMMInputQueue;
    protected LinkedBlockingDeque<AssemblyResultSet> pairHMMResultQueue;
    protected PriorityBlockingQueue<ActiveRegionVariantSet> variantSetQueue;

    protected List<Mutect2AccEngine> dataPrepareEngineList;
    protected List<Mutect2AccEngine> localAssemblyEngineList;
    protected List<Mutect2AccEngine> loglessPairHMMEngineList;
    protected List<Mutect2AccEngine> genotypingEngineList;

    protected Queue<Runnable> dataPrepareRunnableQueue;
    protected Queue<Runnable> localAssemblyRunnableQueue;
    protected Queue<Runnable> loglessPairHMMRunnableQueue;
    protected Queue<Runnable> genotypingRunnableQueue;

    protected Lock dataPrepareRunnableLock;
    protected Lock localAssemblyRunnableLock;
    protected Lock loglessPairHMMRunnableLock;
    protected Lock genotypingRunnableLock;

    protected boolean tooManyEnginesForDataPrepare;
    protected boolean tooManyEnginesForLocalAssembly;
    protected boolean tooManyEnginesForLoglessPairHMM;
    protected boolean tooManyEnginesForGenotyping;

    protected AtomicInteger numActiveRegionCreatedByDataPrepare;
    protected AtomicInteger numActiveRegionAddressedByLocalAssembly;
    protected AtomicInteger numActiveRegionAddressedByPairHMM;
    protected AtomicInteger numActiveRegionAddressedByGenotyping;
    protected AtomicInteger numVariantSet;
    protected AtomicInteger numActiveRegionAssignedToVector;
    protected AtomicInteger numActiveRegionAssignedToLogless;

    protected boolean pairHMMSelectionDone;
    protected boolean vectorLogessPairHMMWorkDone;
    protected boolean runnableIsStarted;

    protected AtomicInteger indexOfCurrentShard;
    protected AtomicInteger numRemainedIntervals;
    protected int[] numActiveRegionOnShards;

    private static final int RESULT_QUEUE_STATUS_LOW = 0;
    private static final int RESULT_QUEUE_STATUS_NORMAL = 1;
    private static final int RESULT_QUEUE_STATUS_HIGH = 2;

    private static final int DATA_PREPARE_RUNNABLE = 0;
    private static final int LOCAL_ASSEMBLY_RUNNABLE = 1;
    private static final int LOGLESS_PAIRHMM_RUNNABLE = 2;
    private static final int GENOTYPING_RUNNABLE = 3;

    private static final int NUM_ENGINE_TYPE = 4;

    private final List<MultiIntervalLocalReadShard> readShards;
    private final CommandLineProgram instance;
    private final CountingReadFilter countedFilter;
    private final ReferenceInputArgumentCollection referenceArguments;
    private final ReadInputArgumentCollection readArguments;
    private final M2ArgumentCollection MTAC;
    private final int cloudPrefetchBuffer;
    private final int cloudIndexPrefetchBuffer;
    private final int maxReadsPerAlignmentStart;
    private final boolean isOMPSupported;
    private final Logger logger;
    private final String tumorSample;
    private final String normalSample;
    private final int minRegionSize;
    private final int maxRegionSize;
    private final int assemblyRegionPadding;
    private final int maxProbPropagationDistance;
    private final double activeProbThreshold;
    private final boolean includeReadsWithDeletionsInIsActivePileups;
    private final M2NumThreads m2threads;
    private final int threadAdjustInterval;
    private final int maxJobQueueLength;
    private final SAMFileHeader header;
    private final SampleList samplesList;
    private final Collection<Annotation> annotations;
    private final boolean lenientVCFProcessing;
    private final VCFHeader vcfHeader;

    private long startTimeTagNano;
    private long startTimeIntervalNano = 300000000000L;
    private long threadAdjustIntervalMillis;
    private long actualThreadAdjustIntervalMillis;
    private boolean isCloseToEnd;
    private int lowerLimitOfNumElementInQueue;
    private int upperLimitOfNumElementInQueue;
    private int halfLimitofNumElementInQueue;
    private int COMPUTATIONAL_INTENSIVE_LIMIT;
    private int MIN_COMPUTATIONAL_INTENSIVE_LIMIT;
    private int ADJUST_STEP;
    private int numPreviousRemainedThreadsForDataPrepare;
    private int numPreviousRemainedThreadsForLocalAssembly;
    private int numPreviousRemainedThreadsForLoglessPairHMM;

    private ReallocateDataSet[] reallocateDataSets;

    private String[] engineNames;

    MonitorRunnable (
            final List<MultiIntervalLocalReadShard> readShards,
            final CommandLineProgram instance,
            final CountingReadFilter countedFilter,
            final ReferenceInputArgumentCollection referenceArguments,
            final ReadInputArgumentCollection readArguments,
            final M2ArgumentCollection MTAC,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer,
            final int maxReadsPerAlignmentStart,
            final boolean isOMPSupported,
            final Logger logger,
            final String tumorSample,
            final String normalSample,
            final int minRegionSize,
            final int maxRegionSize,
            final int assemblyRegionPadding,
            final double activeProbThreshold,
            final int maxProbPropagationDistance,
            final boolean includeReadsWithDeletionsInIsActivePileups,
            final M2NumThreads m2threads,
            final int threadAdjustInterval,
            final int maxJobQueueLength,
            final SAMFileHeader header,
            final SampleList samplesList,
            final Collection<Annotation> annotations,
            final boolean lenientVCFProcessing,
            final VCFHeader vcfHeader
            ) {
        this.readShards = readShards;
        this.instance = instance;
        this.countedFilter = countedFilter;
        this.referenceArguments = referenceArguments;
        this.readArguments = readArguments;
        this.MTAC = MTAC;
        this.cloudPrefetchBuffer = cloudPrefetchBuffer;
        this.cloudIndexPrefetchBuffer = cloudIndexPrefetchBuffer;
        this.maxReadsPerAlignmentStart = maxReadsPerAlignmentStart;
        this.isOMPSupported = isOMPSupported;
        this.logger = logger;
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
        this.minRegionSize = minRegionSize;
        this.maxRegionSize = maxRegionSize;
        this.assemblyRegionPadding = assemblyRegionPadding;
        this.maxProbPropagationDistance = maxProbPropagationDistance;
        this.activeProbThreshold = activeProbThreshold;
        this.includeReadsWithDeletionsInIsActivePileups = includeReadsWithDeletionsInIsActivePileups;
        this.m2threads = m2threads;
        this.threadAdjustInterval = threadAdjustInterval;
        this.maxJobQueueLength = maxJobQueueLength;
        this.header = header;
        this.samplesList = samplesList;
        this.annotations = annotations;
        this.lenientVCFProcessing = lenientVCFProcessing;
        this.vcfHeader = vcfHeader;

        numActiveRegionOnShards = new int[readShards.size()];
        for (int i=0; i<readShards.size(); ++i) {
            numActiveRegionOnShards[i] = -1;
        }

        numRemainedIntervals = new AtomicInteger(readShards.size());

        vectorLogessPairHMMWorkDone = false;
        runnableIsStarted = false;
        indexOfCurrentShard = new AtomicInteger(0);

        initializeJobQueues();

        startTimeTagNano = System.nanoTime();
        threadAdjustIntervalMillis = threadAdjustInterval * 1000;
        isCloseToEnd = false;
        upperLimitOfNumElementInQueue = maxJobQueueLength * 9 / 10;
        lowerLimitOfNumElementInQueue = maxJobQueueLength / 10;
        halfLimitofNumElementInQueue = maxJobQueueLength / 2;

        COMPUTATIONAL_INTENSIVE_LIMIT = 800;
        MIN_COMPUTATIONAL_INTENSIVE_LIMIT = 200;
        ADJUST_STEP = 200;

        numPreviousRemainedThreadsForDataPrepare = m2threads.numThreadsForDataPrepare;
        numPreviousRemainedThreadsForLocalAssembly = m2threads.numThreadsForLocalAssembly;
        numPreviousRemainedThreadsForLoglessPairHMM = m2threads.numThreadsForLoglessPairHMM;

        reallocateDataSets = new ReallocateDataSet[NUM_ENGINE_TYPE];

        reallocateDataSets[DATA_PREPARE_RUNNABLE] = new ReallocateDataSet();
        reallocateDataSets[DATA_PREPARE_RUNNABLE].engineList = dataPrepareEngineList;
        reallocateDataSets[DATA_PREPARE_RUNNABLE].runnableQueue = dataPrepareRunnableQueue;
        reallocateDataSets[DATA_PREPARE_RUNNABLE].runnableLock = dataPrepareRunnableLock;

        reallocateDataSets[LOCAL_ASSEMBLY_RUNNABLE] = new ReallocateDataSet();
        reallocateDataSets[LOCAL_ASSEMBLY_RUNNABLE].engineList = localAssemblyEngineList;
        reallocateDataSets[LOCAL_ASSEMBLY_RUNNABLE].runnableQueue = localAssemblyRunnableQueue;
        reallocateDataSets[LOCAL_ASSEMBLY_RUNNABLE].runnableLock = localAssemblyRunnableLock;

        reallocateDataSets[LOGLESS_PAIRHMM_RUNNABLE] = new ReallocateDataSet();
        reallocateDataSets[LOGLESS_PAIRHMM_RUNNABLE].engineList = loglessPairHMMEngineList;
        reallocateDataSets[LOGLESS_PAIRHMM_RUNNABLE].runnableQueue = loglessPairHMMRunnableQueue;
        reallocateDataSets[LOGLESS_PAIRHMM_RUNNABLE].runnableLock = loglessPairHMMRunnableLock;

        reallocateDataSets[GENOTYPING_RUNNABLE] = new ReallocateDataSet();
        reallocateDataSets[GENOTYPING_RUNNABLE].engineList = genotypingEngineList;
        reallocateDataSets[GENOTYPING_RUNNABLE].runnableQueue = genotypingRunnableQueue;
        reallocateDataSets[GENOTYPING_RUNNABLE].runnableLock =  genotypingRunnableLock;

        engineNames = new String[NUM_ENGINE_TYPE];
        engineNames[DATA_PREPARE_RUNNABLE] = "DATA_PREPARE";
        engineNames[LOCAL_ASSEMBLY_RUNNABLE] = "LOCAL_ASSEMBLY";
        engineNames[LOGLESS_PAIRHMM_RUNNABLE] = "LOGLESS_PAIRHMM";
        engineNames[GENOTYPING_RUNNABLE] = "GENOTYPING";
    }

    public void run () {
        // waiting for initial threads started
        while (numInitializedThreadsForDataPrepare.get() < m2threads.numThreadsForDataPrepare
                || numInitializedThreadsForLocalAssembly.get() < m2threads.numThreadsForLocalAssembly
                || numInitializedThreadsForLoglessPairHMM.get() < m2threads.numThreadsForLoglessPairHMM
                || numInitializedThreadsForGenotyping.get() < m2threads.numThreadsForGenotyping) {
            try {
                Thread.sleep(1);
            } catch (InterruptedException e) {

            }
        }

        System.out.println ("Initial Threads Started, Starting Monitoring");

        while (true) {
            if (numRemainedThreadsForDataPrepare.get() <= 0
                    && numRemainedThreadsForLocalAssembly.get() <= 0
                    && numRemainedThreadsForGenotyping.get() <= 0) {
                break;
            }

            if (System.nanoTime()-startTimeTagNano < startTimeIntervalNano) {
                actualThreadAdjustIntervalMillis = 5000;
            } else if (isCloseToEnd) {
                actualThreadAdjustIntervalMillis = threadAdjustIntervalMillis / 2;
            } else {
                actualThreadAdjustIntervalMillis = threadAdjustIntervalMillis;
            }

            try {
                Thread.sleep(actualThreadAdjustIntervalMillis);
            } catch (InterruptedException e) {

            }

            if (isOMPSupported) {
                int vectorPairHMMInputQueueLength = vectorLoglessPairHMMInputQueue.size();
                int loglessPairHMMInputQueueLength = loglessPairhMMInputQueue.size();
                if (vectorPairHMMInputQueueLength>upperLimitOfNumElementInQueue && loglessPairHMMInputQueueLength<=upperLimitOfNumElementInQueue) {
                    System.out.println("COMPUTATIONAL_INTENSIVE_LIMIT change from " + COMPUTATIONAL_INTENSIVE_LIMIT + " to " + (COMPUTATIONAL_INTENSIVE_LIMIT + ADJUST_STEP));
                    COMPUTATIONAL_INTENSIVE_LIMIT += ADJUST_STEP;
                } else if (vectorPairHMMInputQueueLength<lowerLimitOfNumElementInQueue && loglessPairHMMInputQueueLength>lowerLimitOfNumElementInQueue
                        || vectorPairHMMInputQueueLength<=upperLimitOfNumElementInQueue && loglessPairHMMInputQueueLength>upperLimitOfNumElementInQueue) {
                    if (COMPUTATIONAL_INTENSIVE_LIMIT > MIN_COMPUTATIONAL_INTENSIVE_LIMIT) {
                        System.out.println("COMPUTATIONAL_INTENSIVE_LIMIT change from " + COMPUTATIONAL_INTENSIVE_LIMIT + " to " + (COMPUTATIONAL_INTENSIVE_LIMIT - ADJUST_STEP));
                        COMPUTATIONAL_INTENSIVE_LIMIT -= ADJUST_STEP;
                    }
                }
            }

            System.out.println ("numPreviousRemainedThreadsForDataPrepare: " + numPreviousRemainedThreadsForDataPrepare + "\tnumRemainedThreadsForDataPrepare: " + numRemainedThreadsForDataPrepare.get());
            while (numPreviousRemainedThreadsForDataPrepare > 0
                    && numRemainedThreadsForDataPrepare.get() < numPreviousRemainedThreadsForDataPrepare) {
                // find one available engine, which has finished its job
                if (!isCloseToEnd) {
                    isCloseToEnd = true;
                }

                Mutect2AccEngine availableEngine = null;
                for (final Mutect2AccEngine engine : dataPrepareEngineList) {
                    if (engine.future.isDone()) {
                        availableEngine = engine;
                        break;
                    }
                }
                if (availableEngine == null) {
                    throw new IllegalArgumentException("none of calculation engines are available!");
                }
                dataPrepareEngineList.remove(availableEngine);
                availableEngine.runnable = null;
                availableEngine.future = null;
                --numPreviousRemainedThreadsForDataPrepare;

                // reallocate this engine to LocalAssembly or LoglessPairHMM or Genotyping
                Runnable runnable;
                int pairHMMResultQueueLength = pairHMMResultQueue.size();
                if (numRemainedThreadsForLoglessPairHMM.get() > 0
                        && pairHMMResultQueueLength <= lowerLimitOfNumElementInQueue) {
                    startLoglessPairHMMEngine(availableEngine);
                    ++numPreviousRemainedThreadsForLoglessPairHMM;
                    System.out.println ("Reallocate engine from " + engineNames[DATA_PREPARE_RUNNABLE] + " to " + engineNames[LOGLESS_PAIRHMM_RUNNABLE] + " in step 1");
                } else {
                    startGenotypingEngine(availableEngine);
                    System.out.println ("Reallocate engine from " + engineNames[DATA_PREPARE_RUNNABLE] + " to " + engineNames[GENOTYPING_RUNNABLE] + " in step 1");
                }
            }

            System.out.println ("numPreviousRemainedThreadsForLocalAssembly: " + numPreviousRemainedThreadsForLocalAssembly + "\tnumRemainedThreadsForLocalAssembly: " + numRemainedThreadsForLocalAssembly.get());
            while (numPreviousRemainedThreadsForLocalAssembly > 0
                    && numRemainedThreadsForLocalAssembly.get() < numPreviousRemainedThreadsForLocalAssembly) {
                // Find one available engine, which has finish its job
                if (!isCloseToEnd) {
                    isCloseToEnd = true;
                }

                Mutect2AccEngine availableEngine = null;
                for (final Mutect2AccEngine engine : localAssemblyEngineList) {
                    if (engine.future.isDone()) {
                        availableEngine = engine;
                        break;
                    }
                }
                if (availableEngine == null) {
                    throw new IllegalArgumentException("none of calculaton engines are available");
                }
                localAssemblyEngineList.remove(availableEngine);
                availableEngine.runnable = null;
                availableEngine.future = null;
                --numPreviousRemainedThreadsForLocalAssembly;

                // Reallocate this engine to Logless or Genotyping
                Runnable runnable;
                int pairHMMResultQueueLength = pairHMMResultQueue.size();
                if (numRemainedThreadsForLoglessPairHMM.get() > 0
                        && pairHMMResultQueueLength <= lowerLimitOfNumElementInQueue) {
                    startLoglessPairHMMEngine(availableEngine);
                    ++numPreviousRemainedThreadsForLoglessPairHMM;
                    System.out.println ("Reallocate engine from " + engineNames[LOCAL_ASSEMBLY_RUNNABLE] + " to " + engineNames[LOGLESS_PAIRHMM_RUNNABLE] + " in step 2");
                } else {
                    startGenotypingEngine(availableEngine);
                    System.out.println ("Reallocate engine from " + engineNames[LOCAL_ASSEMBLY_RUNNABLE] + " to " + engineNames[GENOTYPING_RUNNABLE] + " in step 2");
                }
            }

            System.out.println ("numPreviousRemainedThreadsForLoglessPairHMM: " + numPreviousRemainedThreadsForLoglessPairHMM + "\tnumRemainedThreadsForLoglessPairHMM: " + numRemainedThreadsForLoglessPairHMM.get());
            while (numPreviousRemainedThreadsForLoglessPairHMM > 0
                    && numRemainedThreadsForLoglessPairHMM.get() < numPreviousRemainedThreadsForLoglessPairHMM) {
                // Find one available engine, which has finish its job
                if (!isCloseToEnd) {
                    isCloseToEnd = true;
                }

                Mutect2AccEngine availableEngine = null;
                for (final Mutect2AccEngine engine : loglessPairHMMEngineList) {
                    if (engine.future.isDone()) {
                        availableEngine = engine;
                        break;
                    }
                }
                if (availableEngine == null) {
                    throw new IllegalArgumentException("none of calculation engines are available");
                }
                loglessPairHMMEngineList.remove(availableEngine);
                availableEngine.runnable = null;
                availableEngine.future = null;
                --numPreviousRemainedThreadsForLoglessPairHMM;

                // Reallocate this engine to Genotyping
                startGenotypingEngine(availableEngine);
                System.out.println ("Reallocate engine from " + engineNames[LOGLESS_PAIRHMM_RUNNABLE] + " to " + engineNames[GENOTYPING_RUNNABLE] + " in step 3");
            }

            reallocateEngines ();
        }

        // Data Prepare Engines Shutdown
        System.out.println ("Shutdown data prepare engines");
        shutdownEngines (dataPrepareEngineList);

        // Local Assembly Engines Shutdown
        System.out.println ("Shutdown local assembly engines");
        shutdownEngines (localAssemblyEngineList);

        // Logless PairHMM Engines Shutdown
        System.out.println ("Shutdown logless pairHMM engines");
        shutdownEngines (loglessPairHMMEngineList);

        // Genotyping Engines Shutdown
        System.out.println ("Shutdown genotyping engines");
        shutdownEngines (genotypingEngineList);
    }

    private void initializeJobQueues () {
        numRemainedThreadsForDataPrepare = new AtomicInteger(0);
        numRemainedThreadsForLocalAssembly = new AtomicInteger(0);
        numRemainedThreadsForLoglessPairHMM = new AtomicInteger(0);
        numRemainedThreadsForGenotyping = new AtomicInteger(0);

        numInitializedThreadsForDataPrepare = new AtomicInteger(0);
        numInitializedThreadsForLocalAssembly = new AtomicInteger(0);
        numInitializedThreadsForLoglessPairHMM = new AtomicInteger(0);
        numInitializedThreadsForGenotyping = new AtomicInteger(0);

        dataPrepareRunnableLock = new ReentrantLock();
        localAssemblyRunnableLock = new ReentrantLock();
        loglessPairHMMRunnableLock = new ReentrantLock();
        genotypingRunnableLock = new ReentrantLock();

        tooManyEnginesForDataPrepare = false;
        tooManyEnginesForLocalAssembly = false;
        tooManyEnginesForLoglessPairHMM = false;
        tooManyEnginesForGenotyping = false;

        numActiveRegionCreatedByDataPrepare = new AtomicInteger(0);
        numActiveRegionAddressedByLocalAssembly = new AtomicInteger(0);
        numActiveRegionAddressedByPairHMM = new AtomicInteger(0);
        numActiveRegionAddressedByGenotyping = new AtomicInteger(0);
        numVariantSet = new AtomicInteger(0);

        numActiveRegionAssignedToVector = new AtomicInteger(0);
        numActiveRegionAssignedToLogless = new AtomicInteger(0);

        activeRegionQueue = new LinkedBlockingDeque<>(Mutect2Acc.MAX_QUEUE_LENGTH);
        localAssemblyResultQueue = new LinkedBlockingDeque<>(Mutect2Acc.MAX_QUEUE_LENGTH);
        pairHMMResultQueue = new LinkedBlockingDeque<>(Mutect2Acc.MAX_QUEUE_LENGTH);
        if (isOMPSupported) {
            vectorLoglessPairHMMInputQueue = new LinkedBlockingDeque<>(Mutect2Acc.MAX_QUEUE_LENGTH);
            loglessPairhMMInputQueue = new LinkedBlockingDeque<>(Mutect2Acc.MAX_QUEUE_LENGTH);
        }
        variantSetQueue = new PriorityBlockingQueue<>();

        dataPrepareEngineList = new LinkedList<>();
        localAssemblyEngineList = new LinkedList<>();
        loglessPairHMMEngineList = new LinkedList<>();
        genotypingEngineList = new LinkedList<>();

        dataPrepareRunnableQueue = new LinkedList<>();
        localAssemblyRunnableQueue = new LinkedList<>();
        loglessPairHMMRunnableQueue = new LinkedList<>();
        genotypingRunnableQueue = new LinkedList<>();
    }

    private void reallocateEngines () {
        int activeRegionQueueLength = activeRegionQueue.size();
        int localAssemblyQueueLength = localAssemblyResultQueue.size();
        int pairHMMQueueLength = pairHMMResultQueue.size();
        int variantQueueLength = variantSetQueue.size();

        int activeRegionQueueStatus = getQueueStatus(activeRegionQueueLength);
        int localAssemblyQueueStatus = getQueueStatus(localAssemblyQueueLength);
        int pairHMMQueueStatus = getQueueStatus(pairHMMQueueLength);

        int vectorPairHMMInputQueueLength = -1;
        int loglessPairHMMInputQueueLength = -1;
        boolean vectorPairHMMInputQueueMoreThanHalf = false;
        boolean loglessPairHMMInputQueueMoreThanHalf = false;

        if (isOMPSupported) {
            vectorPairHMMInputQueueLength = vectorLoglessPairHMMInputQueue.size();
            loglessPairHMMInputQueueLength = loglessPairhMMInputQueue.size();
            vectorPairHMMInputQueueMoreThanHalf = vectorPairHMMInputQueueLength > halfLimitofNumElementInQueue;
            loglessPairHMMInputQueueMoreThanHalf = loglessPairHMMInputQueueLength > halfLimitofNumElementInQueue;
        }

        System.out.printf ("Accumulated Time: %fm\n", (double)(System.nanoTime()-startTimeTagNano) / 60000000000.0);
        System.out.printf ("Remained interval cnt: %d\tData cnt: %d\tAssembly result cnt: %d\tPairHMM result cnt: %d\tVariant set cnt: %d\n",
                numRemainedIntervals.get(), activeRegionQueueLength, localAssemblyQueueLength, pairHMMQueueLength, variantQueueLength);
        if (isOMPSupported) {
            System.out.printf ("Vector pairHMM input cnt: %d\tLogless pairHMM input cnt:%d\n",
                    vectorPairHMMInputQueueLength, loglessPairHMMInputQueueLength);
            System.out.printf ("Regions allocated to vector pairHMM: %d\tRegions allocated to logless pairHMM: %d\n",
                    numActiveRegionAssignedToVector.get(), numActiveRegionAssignedToLogless.get());
        }
        System.out.printf ("nth for DP: %d\tnth for LA: %d\tnth for LP: %d\tnth for GT: %d\n",
                numRemainedThreadsForDataPrepare.get(), numRemainedThreadsForLocalAssembly.get(),
                numRemainedThreadsForLoglessPairHMM.get(), numRemainedThreadsForGenotyping.get());
        System.out.printf ("DP AR: %d\tLA AR: %d\tPH: AR: %d\tGT AR: %d\n\n",
                numActiveRegionCreatedByDataPrepare.get(), numActiveRegionAddressedByLocalAssembly.get(),
                numActiveRegionAddressedByPairHMM.get(), numActiveRegionAddressedByGenotyping.get());
        System.out.flush();

        switch (activeRegionQueueStatus) {
            case RESULT_QUEUE_STATUS_LOW:
                switch (localAssemblyQueueStatus) {
                    case RESULT_QUEUE_STATUS_LOW:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    if (vectorPairHMMInputQueueMoreThanHalf && loglessPairHMMInputQueueMoreThanHalf) {
                                       boolean isSucceeded = reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                       if (!isSucceeded) {
                                           isSucceeded = reallocateEngineCore(GENOTYPING_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                           if (!isSucceeded) {
                                               reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                           }
                                       }
                                    } else {
                                        boolean isSucceeded = reallocateEngineCore(GENOTYPING_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                        if (!isSucceeded) {
                                            isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                            if (!isSucceeded) {
                                                reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                            }
                                        }
                                    }
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (!isCloseToEnd) {
                                    if (vectorPairHMMInputQueueMoreThanHalf && loglessPairHMMInputQueueMoreThanHalf) {
                                        reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                    } else {
                                        boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                        if (!isSucceeded) {
                                            reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                        }
                                    }
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                } else {
                                    boolean isSucceeded = reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    }
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_NORMAL:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_HIGH:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, DATA_PREPARE_RUNNABLE, false);
                                }
                                break;
                        }
                        break;
                }
                break;
            case RESULT_QUEUE_STATUS_NORMAL:
                switch (localAssemblyQueueStatus) {
                    case RESULT_QUEUE_STATUS_LOW:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    boolean isSucceeded = reallocateEngineCore(GENOTYPING_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(GENOTYPING_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_NORMAL:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(GENOTYPING_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_HIGH:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                break;
                        }
                        break;
                }
                break;
            case RESULT_QUEUE_STATUS_HIGH:
                switch (localAssemblyQueueStatus) {
                    case RESULT_QUEUE_STATUS_LOW:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    boolean isSucceeded = reallocateEngineCore(GENOTYPING_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                } else {
                                    boolean isSucceeded = reallocateEngineCore(GENOTYPING_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                } else {
                                    boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(GENOTYPING_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                        reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                } else {
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_NORMAL:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                        reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                } else {
                                    boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    }
                                }
                                break;
                        }
                        break;
                    case RESULT_QUEUE_STATUS_HIGH:
                        switch (pairHMMQueueStatus) {
                            case RESULT_QUEUE_STATUS_LOW:
                                if (!isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, true);
                                    reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, true);
                                } else {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_NORMAL:
                                if (isCloseToEnd) {
                                    reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, LOGLESS_PAIRHMM_RUNNABLE, false);
                                }
                                break;
                            case RESULT_QUEUE_STATUS_HIGH:
                                if (!isCloseToEnd) {
                                    boolean isSucceeded = reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                        reallocateEngineCore(DATA_PREPARE_RUNNABLE, LOCAL_ASSEMBLY_RUNNABLE, false);
                                    }
                                } else {
                                    boolean isSucceeded = reallocateEngineCore(LOCAL_ASSEMBLY_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    if (!isSucceeded) {
                                        reallocateEngineCore(LOGLESS_PAIRHMM_RUNNABLE, GENOTYPING_RUNNABLE, false);
                                    }
                                }
                                break;
                        }
                        break;
                }
                break;
        }
    }

    private boolean reallocateEngineCore (final int srcRunnableType, final int dstRunnableType, final boolean force) {
        List<Mutect2AccEngine> srcEngineList = reallocateDataSets[srcRunnableType].engineList;
        List<Mutect2AccEngine> dstEngineList = reallocateDataSets[dstRunnableType].engineList;
        Queue<Runnable> srcRunnableQueue = reallocateDataSets[srcRunnableType].runnableQueue;
        Queue<Runnable> dstRunnableQueue = reallocateDataSets[dstRunnableType].runnableQueue;
        Lock srcRunnableLock = reallocateDataSets[srcRunnableType].runnableLock;

        switch (dstRunnableType) {
            case DATA_PREPARE_RUNNABLE:
                if (!force && numRemainedThreadsForDataPrepare.get()>=m2threads.maxNumThreadsForDataPrepare) {
                    return false;
                }
                break;
            case LOCAL_ASSEMBLY_RUNNABLE:
                if (!force && numRemainedThreadsForLocalAssembly.get()>=m2threads.maxNumThreadsForLocalAssembly) {
                    return false;
                }
                break;
            case LOGLESS_PAIRHMM_RUNNABLE:
                if (!force && numRemainedThreadsForLoglessPairHMM.get()>=m2threads.maxNumThreadsForLoglessPairHMM) {
                    return false;
                }
                break;
            case GENOTYPING_RUNNABLE:
                if (!force && numRemainedThreadsForGenotyping.get()>=m2threads.maxNumThreadsForGenotyping) {
                    return false;
                }
                break;
        }

        srcRunnableLock.lock();
        switch (srcRunnableType) {
            case DATA_PREPARE_RUNNABLE:
                if ((!force && numRemainedThreadsForDataPrepare.get()<=m2threads.minNumThreadsForDataPrepare)
                        || (force && numRemainedThreadsForDataPrepare.get()<=1)) {
                    srcRunnableLock.unlock();
                    return false;
                }
                tooManyEnginesForDataPrepare = true;
                --numPreviousRemainedThreadsForDataPrepare;
                break;
            case LOCAL_ASSEMBLY_RUNNABLE:
                if ((!force && numRemainedThreadsForLocalAssembly.get()<=m2threads.minNumThreadsForLocalAssembly)
                        || (force && numRemainedThreadsForLocalAssembly.get()<=1)) {
                    srcRunnableLock.unlock();
                    return false;
                }
                tooManyEnginesForLocalAssembly = true;
                --numPreviousRemainedThreadsForLocalAssembly;
                break;
            case LOGLESS_PAIRHMM_RUNNABLE:
                if ((!force && numRemainedThreadsForLoglessPairHMM.get()<=m2threads.minNumThreadsForLoglessPairHMM)
                        || (force && numRemainedThreadsForLoglessPairHMM.get()<=1)) {
                    srcRunnableLock.unlock();
                    return false;
                }
                tooManyEnginesForLoglessPairHMM = true;
                --numPreviousRemainedThreadsForLoglessPairHMM;
                break;
            case GENOTYPING_RUNNABLE:
                if ((!force && numRemainedThreadsForGenotyping.get()<=m2threads.minNumThreadsForGenotyping)
                        || (force && numRemainedThreadsForGenotyping.get()<=1)) {
                    srcRunnableLock.unlock();
                    return false;
                }
                tooManyEnginesForGenotyping = true;
                break;
        }
        srcRunnableLock.unlock();

        Mutect2AccEngine availableEngine = null;
        while (true) {
            for (final Mutect2AccEngine engine : srcEngineList) {
                if (engine.future.isDone()) {
                    availableEngine = engine;
                    break;
                }
            }

            if (availableEngine != null) {
                break;
            }

            try {
                Thread.sleep(1);
            } catch (InterruptedException e) {

            }
        }
        srcEngineList.remove(availableEngine);
        availableEngine.future = null;
        srcRunnableQueue.add (availableEngine.runnable);

        switch (dstRunnableType) {
            case DATA_PREPARE_RUNNABLE:
                startDataPrepareEngine(availableEngine);
                ++numPreviousRemainedThreadsForDataPrepare;
                break;
            case LOCAL_ASSEMBLY_RUNNABLE:
                startLocalAssemblyEngine(availableEngine);
                ++numPreviousRemainedThreadsForLocalAssembly;
                break;
            case LOGLESS_PAIRHMM_RUNNABLE:
                startLoglessPairHMMEngine(availableEngine);
                ++numPreviousRemainedThreadsForLoglessPairHMM;
                break;
            case GENOTYPING_RUNNABLE:
                startGenotypingEngine(availableEngine);
                break;
        }

        System.out.println ("Reallocate engine from " + engineNames[srcRunnableType] + " to " + engineNames[dstRunnableType] + " in step 4");

        return true;
    }

    private void startDataPrepareEngine (final Mutect2AccEngine engine) {
       Runnable runnable = dataPrepareRunnableQueue.poll();
       if (runnable == null) {
           runnable = new DataPrepareRunnable(
                   readShards,
                   instance,
                   this,
                   countedFilter,
                   referenceArguments,
                   readArguments,
                   MTAC,
                   cloudPrefetchBuffer,
                   cloudIndexPrefetchBuffer,
                   maxReadsPerAlignmentStart,
                   logger,
                   tumorSample,
                   normalSample,
                   minRegionSize,
                   maxRegionSize,
                   assemblyRegionPadding,
                   activeProbThreshold,
                   maxProbPropagationDistance,
                   includeReadsWithDeletionsInIsActivePileups
           );
       }
       startEngineCore(engine, runnable, dataPrepareEngineList);
    }

    private void startLocalAssemblyEngine (final Mutect2AccEngine engine) {
        Runnable runnable = localAssemblyRunnableQueue.poll();
        if (runnable == null) {
            runnable = new LocalAssemblyRunnable(
                    this,
                    MTAC,
                    referenceArguments,
                    instance,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer,
                    header,
                    samplesList,
                    logger
            );
        }
        startEngineCore(engine, runnable, localAssemblyEngineList);
    }

    private void startLoglessPairHMMEngine (final Mutect2AccEngine engine) {
        Runnable runnable = loglessPairHMMRunnableQueue.poll();
        if (runnable == null) {
            runnable = new LoglessPairHMMRunnable(this, isOMPSupported, MTAC, samplesList);
        }
        startEngineCore(engine, runnable, loglessPairHMMEngineList);
    }

    private void startGenotypingEngine (final Mutect2AccEngine engine) {
        Runnable runnable = genotypingRunnableQueue.poll();
        if (runnable == null) {
            runnable = new GenotypingRunnable(
                    this,
                    instance,
                    MTAC,
                    referenceArguments,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer,
                    samplesList,
                    tumorSample,
                    normalSample,
                    annotations,
                    header,
                    vcfHeader,
                    lenientVCFProcessing
            );
        }
        startEngineCore(engine, runnable, genotypingEngineList);
    }

    private void startEngineCore (final Mutect2AccEngine engine, final Runnable runnable, List<Mutect2AccEngine> engineList) {
        runnableIsStarted = false;
        engine.runnable = runnable;
        engine.future = engine.threadPool.submit(runnable);
        engineList.add (engine);
        while (!runnableIsStarted) {
            try {
                Thread.sleep(1);
            } catch (InterruptedException e) {

            }
        }
    }

    private int getQueueStatus (final int queueLength) {
        if (queueLength < lowerLimitOfNumElementInQueue) {
            return RESULT_QUEUE_STATUS_LOW;
        }

        if (queueLength < upperLimitOfNumElementInQueue) {
            return RESULT_QUEUE_STATUS_NORMAL;
        }

        return RESULT_QUEUE_STATUS_HIGH;
    }

    private void shutdownEngines (final List<Mutect2AccEngine> engines) {
       for (final Mutect2AccEngine engine : engines) {
           final ExecutorService threadPool = engine.threadPool;
           threadPool.shutdown();
           while (true) {
               if (threadPool.isTerminated()) {
                   break;
               }

               try {
                   Thread.sleep(1);
               } catch (InterruptedException e) {

               }
           }
       }
    }

    private static class ReallocateDataSet {
        public List<Mutect2AccEngine> engineList;
        public Queue<Runnable> runnableQueue;
        public Lock runnableLock;
    }
}
