package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pairhmm.VectorLoglessPairHMM;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

@CommandLineProgramProperties(
        summary = "Fast call somatic SNVs and indels via local assembly of haplotypes",
        oneLineSummary = "Fast call somatic SNVs and indels via local assembly of haplotypes",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)

@DocumentedFeature
public final class Mutect2Acc extends AssemblyRegionWalker {

    @ArgumentCollection
    private M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    private File outputVCF;

    @Argument(fullName = "num-threads", doc = "How many threads Mutect2 use (default: max-cpu-threads/2)", optional = true)
    private int totalNumThreads = -1;

    @Argument(fullName = "thread-adjust-interval", doc = "How many seconds / per thread adjust (default: 30)", optional = true)
    private int threadAdjustInterval = 30;

    @Argument(fullName = "max-job-queue-length", doc = "Max length of the job queues, must >= 256 (default: 1024)", optional = true)
    static protected int MAX_QUEUE_LENGTH = 4096;

    @Argument(fullName = "use-fpga", doc = "Whether to use FPGA or not", optional = true)
    public static boolean useFPGA = false;

    @Argument(fullName = "native-pair-hmm-threads", doc="How many threads should a native pairHMM implementation use (default: max-cpu-threads/2)", optional = true)
    public static int pairHmmNativeThreads = -1;

    private SAMFileHeader header;
    private SAMSequenceDictionary sequenceDictionary;
    private VariantContextWriter vcfWriter;
    private CountingReadFilter countedFilter;
    private VCFHeader vcfHeader;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private Collection<Annotation> annotations;

    private SampleList samplesList;
    private String tumorSample;
    private String normalSample;

    private M2NumThreads m2threads;

    private boolean isOMPSupported;

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer().andThen(Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), !MTAC.dontClipITRArtifacts));
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Mutect2Engine.getStandardMutect2AnnotationGroups();
    }

    @Override
    protected ReadsDownsampler createDownsampler() {
        return new MutectDownsampler(maxReadsPerAlignmentStart, MTAC.maxSuspiciousReadsPerAlignmentStart, MTAC.downsamplingStride);
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        return null;
    }

    @Override
    public void onTraversalStart() {
        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? userIntervals : loadEffectiveRegionForReference(referenceArguments.getReferenceFileName());
        readShards = makeReadShards (intervals);

        annotations = makeVariantAnnotations();
        header = reads.getHeader();
        sequenceDictionary = reads.getSequenceDictionary();
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(header)));
        tumorSample = decodeSampleNameIfNecessary(MTAC.tumorSample);
        normalSample = MTAC.normalSample==null ? null : decodeSampleNameIfNecessary(MTAC.normalSample);

        vcfWriter = createVCFWriter(outputVCF);
        writeVCFHeader(getDefaultToolVCFHeaderLines());

        if (MAX_QUEUE_LENGTH < 256) {
            MAX_QUEUE_LENGTH = 256;
        }

        countedFilter = makeReadFilter();

        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs);
        isOMPSupported = checkOMP();

        initializeThreads ();
    }

    @Override
    public Object onTraversalSuccess () {
        return "SUCCESS";
    }

    @Override
    public void apply (final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext) {

    }

    @Override
    public final void traverse() {
        // Monitor Start
        MonitorRunnable monitor = new MonitorRunnable(
                readShards,
                this,
                countedFilter,
                referenceArguments,
                readArguments,
                MTAC,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer,
                maxReadsPerAlignmentStart,
                isOMPSupported,
                logger,
                tumorSample,
                normalSample,
                minAssemblyRegionSize,
                maxAssemblyRegionSize,
                assemblyRegionPadding,
                activeProbThreshold,
                maxProbPropagationDistance,
                includeReadsWithDeletionsInIsActivePileups(),
                m2threads,
                threadAdjustInterval,
                MAX_QUEUE_LENGTH,
                header,
                samplesList,
                annotations,
                lenientVCFProcessing,
                vcfHeader
        );
        ExecutorService monitorThreadPool = Executors.newSingleThreadExecutor();
        monitorThreadPool.execute(monitor);

        // Active Region Preparation Engines Start
        for (int i=0; i<m2threads.numThreadsForDataPrepare; ++i) {
            DataPrepareRunnable runnable = new DataPrepareRunnable(
                    readShards,
                    this,
                    monitor,
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
                    minAssemblyRegionSize,
                    maxAssemblyRegionSize,
                    assemblyRegionPadding,
                    activeProbThreshold,
                    maxProbPropagationDistance,
                    includeReadsWithDeletionsInIsActivePileups()
            );
            ExecutorService threadPool = Executors.newSingleThreadExecutor();
            Future<?> future = threadPool.submit(runnable);
            monitor.dataPrepareEngineList.add (new Mutect2AccEngine(threadPool,runnable,future));
        }

        // Local Assembly Engines Start
        for (int i=0; i<m2threads.numThreadsForLocalAssembly; ++i) {
            LocalAssemblyRunnable runnable = new LocalAssemblyRunnable(
                    monitor,
                    MTAC,
                    referenceArguments,
                    this,
                    cloudPrefetchBuffer,
                    cloudIndexPrefetchBuffer,
                    header,
                    samplesList,
                    logger
            );
            ExecutorService threadPool = Executors.newSingleThreadExecutor();
            Future<?> future = threadPool.submit(runnable);
            monitor.localAssemblyEngineList.add (new Mutect2AccEngine(threadPool,runnable,future));
        }

        // PairHMM Selector Engine Start
        ExecutorService pairHMMThreadPool = null;
        if (isOMPSupported) {
            PairHMMJobSelectorRunnable runnable = new PairHMMJobSelectorRunnable(monitor);
            pairHMMThreadPool = Executors.newSingleThreadExecutor();
            pairHMMThreadPool.execute(runnable);
        }

        // Logless PairHMM Engines Start
        for (int i=0; i<m2threads.numThreadsForLoglessPairHMM; ++i) {
            LoglessPairHMMRunnable runnable = new LoglessPairHMMRunnable(
                    monitor,
                    isOMPSupported,
                    MTAC,
                    samplesList
            );
            ExecutorService threadPool = Executors.newSingleThreadExecutor();
            Future<?> future = threadPool.submit(runnable);
            monitor.loglessPairHMMEngineList.add (new Mutect2AccEngine(threadPool,runnable,future));
        }

        // Genotyping Engines Start
        for (int i=0; i<m2threads.numThreadsForGenotyping; ++i) {
            GenotypingRunnable runnable = new GenotypingRunnable(
                    monitor,
                    this,
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
            ExecutorService threadPool = Executors.newSingleThreadExecutor();
            Future<?> future = threadPool.submit(runnable);
            monitor.genotypingEngineList.add (new Mutect2AccEngine(threadPool,runnable,future));
        }

        // OMP/FPGA PairHMM Calculation
        if (isOMPSupported) {
            while (true) {
                AssemblyResultSet assemblyResult = monitor.vectorLoglessPairHMMInputQueue.poll();
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

                assemblyResult.readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResult, samplesList, assemblyResult.reads);

                monitor.numActiveRegionAddressedByPairHMM.incrementAndGet();

                while (!monitor.pairHMMResultQueue.offer(assemblyResult)) {
                    try {
                        Thread.sleep(1);
                    } catch (InterruptedException e) {

                    }
                }
            }
        }
        monitor.vectorLogessPairHMMWorkDone = true;

        if (pairHMMThreadPool != null) {
            pairHMMThreadPool.shutdown();
            while (!pairHMMThreadPool.isTerminated()) {
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
            }
        }

        monitorThreadPool.shutdown();
        while (!monitorThreadPool.isTerminated()) {
            try {
                Thread.sleep(1);
            } catch (InterruptedException e) {

            }
        }

        // Dump Variants
        long beg = System.nanoTime();
        ActiveRegionVariantSet set;
        while ((set = monitor.variantSetQueue.poll()) != null) {
            for (final VariantContext var : set.vaiants) {
                vcfWriter.dump(var);
            }
        }
        System.out.printf ("Variant Dumping Cost: %fs\n", (double)(System.nanoTime()-beg)/1000000000.0);

        logger.info(countedFilter.getSummaryLine());
    }

    @Override
    public void closeTool() {
        vcfWriter.close();
    }

    private List<SimpleInterval> loadEffectiveRegionForReference(final String referencePath) {
        final File effectiveBedFile = new File(referencePath.concat(".effective.bed"));
        BufferedReader effectiveBedReader = null;
        List<SimpleInterval> intervals = new ArrayList<>();

        try {
            effectiveBedReader = new BufferedReader(new FileReader(effectiveBedFile));
        } catch (FileNotFoundException e) {

        }

        try {
            String line = null;
            while ((line = effectiveBedReader.readLine()) != null) {
                final String[] subStrings = line.split("\t");
                // TODO: length of subStrings should be checked
                intervals.add (new SimpleInterval(subStrings[0],
                        Integer.parseInt(subStrings[1])+1, Integer.parseInt(subStrings[2])));
            }
        } catch (IOException e) {

        }

        return intervals;
    }

    private List<MultiIntervalLocalReadShard> makeReadShards(final List<SimpleInterval> intervals ) {
        final List<MultiIntervalLocalReadShard> shards = new ArrayList<>();
        for (final SimpleInterval interval : intervals) {
            shards.add(new MultiIntervalLocalReadShard(Arrays.asList(interval), assemblyRegionPadding, reads.getSequenceDictionary()));
        }

        return shards;
    }

    private String decodeSampleNameIfNecessary (final String name) {
        return samplesList.asListOfSamples().contains(name) ? name : IOUtils.urlDecode(name);
    }

    private void writeVCFHeader (final Set<VCFHeaderLine> defaultToolHeaderLines) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();

        headerInfo.add(new VCFHeaderLine("Mutect Version", Mutect2Engine.MUTECT_VERSION));
        headerInfo.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "Warning: unfiltered Mutect2 calls. Please run " + FilterMutectCalls.class.getSimpleName() + " to remove false positives."));

        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false);
        headerInfo.addAll(annotatorEngine.getVCFAnnotationDescriptions(false));
        headerInfo.addAll(defaultToolHeaderLines);

        GATKVCFConstants.STANDARD_MUTECT_INFO_FIELDS.stream().map(GATKVCFHeaderLines::getInfoLine).forEach(headerInfo::add);

        VCFStandardHeaderLines.addStandardFormatLines(
                headerInfo,
                true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.ALLELE_FRACTION_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY));
        headerInfo.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY));

        headerInfo.add(new VCFHeaderLine(Mutect2Engine.TUMOR_SAMPLE_KEY_IN_VCF_HEADER, tumorSample));
        if (normalSample != null) {
            headerInfo.add(new VCFHeaderLine(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER, normalSample));
        }

        vcfHeader = new VCFHeader(headerInfo, samplesList.asListOfSamples());
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        vcfWriter.writeHeader(vcfHeader);
    }

    private boolean checkOMP () {
        return likelihoodCalculationEngine.getRunningEngineType() == VectorLoglessPairHMM.Implementation.FPGA
                || likelihoodCalculationEngine.getRunningEngineType() == VectorLoglessPairHMM.Implementation.OMP;
    }

    private void initializeThreads () {
        if (totalNumThreads < 0) {
            if (likelihoodCalculationEngine.getRunningEngineType() == VectorLoglessPairHMM.Implementation.FPGA) {
                totalNumThreads = Runtime.getRuntime().availableProcessors();
            } else if (likelihoodCalculationEngine.getRunningEngineType() == VectorLoglessPairHMM.Implementation.OMP) {
                totalNumThreads = Runtime.getRuntime().availableProcessors() - pairHmmNativeThreads;
            } else {
                totalNumThreads = Runtime.getRuntime().availableProcessors();
            }
        }

        m2threads = new M2NumThreads();

        if (isOMPSupported) {
            m2threads.minNumThreadsForDataPrepare = totalNumThreads * 6 / 24;
        } else {
            m2threads.minNumThreadsForDataPrepare = totalNumThreads * 2 / 24;
        }
        m2threads.minNumThreadsForLocalAssembly = totalNumThreads * 3 / 24;
        m2threads.minNumThreadsForLoglessPairHMM = 1;
        m2threads.minNumThreadsForGenotyping = 1;

        m2threads.maxNumThreadsForDataPrepare = totalNumThreads * 18 / 24;
        m2threads.maxNumThreadsForLocalAssembly = totalNumThreads;
        m2threads.maxNumThreadsForLoglessPairHMM = totalNumThreads;
        m2threads.maxNumThreadsForGenotyping = totalNumThreads;

        if (m2threads.minNumThreadsForDataPrepare <= 0) m2threads.minNumThreadsForDataPrepare = 1;
        if (m2threads.minNumThreadsForLocalAssembly <= 0) m2threads.minNumThreadsForLocalAssembly = 1;
        if (m2threads.maxNumThreadsForDataPrepare <= 0) m2threads.maxNumThreadsForDataPrepare = 1;

        if (isOMPSupported) {
            m2threads.numThreadsForLocalAssembly = totalNumThreads * 12 / 24;
            m2threads.numThreadsForLoglessPairHMM = totalNumThreads * 2 / 24;
            m2threads.numThreadsForGenotyping = totalNumThreads * 2 / 24;
            m2threads.numThreadsForDataPrepare = totalNumThreads - m2threads.numThreadsForLocalAssembly - m2threads.numThreadsForLoglessPairHMM - m2threads.numThreadsForGenotyping;
        } else {
            m2threads.numThreadsForDataPrepare = totalNumThreads * 2 / 24;
            m2threads.numThreadsForLocalAssembly = totalNumThreads * 6 / 24;
            m2threads.numThreadsForGenotyping = totalNumThreads * 2 / 24;
            m2threads.numThreadsForLoglessPairHMM = totalNumThreads - m2threads.numThreadsForDataPrepare - m2threads.numThreadsForLocalAssembly - m2threads.numThreadsForGenotyping;
        }

        if (m2threads.numThreadsForDataPrepare <= 0) m2threads.numThreadsForDataPrepare = 1;
        if (m2threads.numThreadsForLocalAssembly <= 0) m2threads.numThreadsForLocalAssembly = 1;
        if (m2threads.numThreadsForLoglessPairHMM <= 0) m2threads.numThreadsForLoglessPairHMM = 1;
        if (m2threads.numThreadsForGenotyping <= 0) m2threads.numThreadsForGenotyping = 1;

        System.out.println();
        System.out.printf ("Total Number of Threads: %d\n", totalNumThreads);
        System.out.println();

        System.out.printf ("Initial Number of Threads for Data Preparation: %d\n", m2threads.numThreadsForDataPrepare);
        System.out.printf ("Initial Number of Threads for Local Assembly: %d\n", m2threads.numThreadsForLocalAssembly);
        System.out.printf ("Initial Number of Threads for Logless PairHMM: %d\n", m2threads.numThreadsForLoglessPairHMM);
        System.out.printf ("Initial Number of Threads for Genotyping: %d\n", m2threads.numThreadsForGenotyping);
        System.out.println();

        System.out.printf ("Minimum Number of Threads for Data Preparation: %d\n", m2threads.minNumThreadsForDataPrepare);
        System.out.printf ("Minimum Number of Threads for Local Assembly: %d\n", m2threads.minNumThreadsForLocalAssembly);
        System.out.printf ("Minimum Number of Threads for Logless PairHMM: %d\n", m2threads.minNumThreadsForLoglessPairHMM);
        System.out.printf ("Minimum Number of Threads for Genotyping: %d\n", m2threads.minNumThreadsForGenotyping);
        System.out.println();

        System.out.printf ("Maximum Number of Threads for Data Preparation: %d\n", m2threads.maxNumThreadsForDataPrepare);
        System.out.printf ("Maximum Number of Threads for Local Assembly: %d\n", m2threads.maxNumThreadsForLocalAssembly);
        System.out.printf ("Maximum Number of Threads for Logless PairHMM: %d\n", m2threads.maxNumThreadsForLoglessPairHMM);
        System.out.printf ("Maximum Number of Threads for Genotyping: %d\n", m2threads.maxNumThreadsForGenotyping);
        System.out.println();
    }



    static protected ReadsDataSource LocalInitializeReads(
            final ReferenceInputArgumentCollection referenceArguments,
            final ReadInputArgumentCollection readArguments,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer
    ) {
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
        factory = factory.referenceSequence(referenceArguments.getReferencePath());
        return new ReadsDataSource(
                readArguments.getReadPaths(),
                readArguments.getReadIndexPaths(),
                factory,
                cloudPrefetchBuffer,
                (cloudIndexPrefetchBuffer<0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer)
        );
    }

    static protected ReferenceDataSource LocalInitializeReference(
            final ReferenceInputArgumentCollection referenceArguments
    ) {
        return referenceArguments.getReferenceFileName()!=null ? ReferenceDataSource.of(referenceArguments.getReferencePath()) : null;
    }

    static protected FeatureManager LocalInitializeFeatures(
            final CommandLineProgram toolInstance,
            final ReferenceInputArgumentCollection referenceArguments,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer
    ) {
        FeatureManager features = new FeatureManager(
                toolInstance,
                FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer,
                referenceArguments.getReferencePath());

        if (features.isEmpty())
            features = null;

        return features;
    }
}
