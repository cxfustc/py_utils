package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.Collections;
import java.util.List;
import java.util.SortedSet;
import java.util.stream.Collectors;

public class LocalAssemblyRunnable implements Runnable {
    private final MonitorRunnable monitor;
    private final M2ArgumentCollection MTAC;
    private final SAMFileHeader header;
    private final SampleList samplesList;
    private final Logger logger;

    private static final int READ_QUALITY_FILTER_THRESHOLD = 20;
    private static final int MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;

    private FeatureManager features;
    private CachingIndexedFastaSequenceFile referenceReader;
    private ReadThreadingAssembler assemblyEngine;
    private SmithWatermanAligner aligner;
    private AssemblyRegionTrimmer trimmer;

    LocalAssemblyRunnable (
            final MonitorRunnable monitor,
            final M2ArgumentCollection MTAC,
            final ReferenceInputArgumentCollection referenceArguments,
            final CommandLineProgram instance,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer,
            final SAMFileHeader header,
            final SampleList samplesList,
            final Logger logger
    ) {
       this.monitor = monitor;
       this.MTAC = MTAC;
       this.header = header;
       this.samplesList = samplesList;
       this.logger = logger;

       features = Mutect2Acc.LocalInitializeFeatures(instance, referenceArguments, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
       referenceReader = AssemblyBasedCallerUtils.createReferenceReader(Utils.nonNull(referenceArguments.getReferenceFileName()));
       assemblyEngine = AssemblyBasedCallerUtils.createReadThreadingAssembler(MTAC);
       aligner = SmithWatermanAligner.getAligner(MTAC.smithWatermanImplementation);
       trimmer = new AssemblyRegionTrimmer();
       trimmer.initialize(
               MTAC.assemblyRegionTrimmerArgs,
               header.getSequenceDictionary(),
               MTAC.debug,
               MTAC.genotypingOutputMode==GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES,
               false
       );
    }

    public void run () {
       monitor.numInitializedThreadsForLocalAssembly.incrementAndGet();

       monitor.localAssemblyRunnableLock.lock();
       monitor.numRemainedThreadsForLocalAssembly.incrementAndGet();
       monitor.localAssemblyRunnableLock.unlock();

       monitor.runnableIsStarted = true;

       while (true) {
           AssemblyRegion assemblyRegion = monitor.activeRegionQueue.poll();
           if (assemblyRegion == null) {
               if (monitor.numRemainedThreadsForDataPrepare.get() <= 0) {
                   break;
               }

               try {
                   Thread.sleep(1);
               } catch (InterruptedException e) {

               }

               if (monitor.tooManyEnginesForLocalAssembly && monitor.localAssemblyRunnableLock.tryLock()) {
                   if (monitor.tooManyEnginesForLocalAssembly) {
                       monitor.tooManyEnginesForLocalAssembly = false;
                       monitor.localAssemblyRunnableLock.unlock();
                       break;
                   }
                   monitor.localAssemblyRunnableLock.unlock();
               }

               continue;
           }

           AssemblyResultSet assemblyResult = doWork (
                   assemblyRegion,
                   new FeatureContext(features, assemblyRegion.getExtendedSpan())
           );
           monitor.numActiveRegionAddressedByLocalAssembly.incrementAndGet();

           if (assemblyResult == null) {
               continue;
           }

           while (!monitor.localAssemblyResultQueue.offer(assemblyResult)) {
               try {
                   Thread.sleep(1);
               } catch (InterruptedException e) {

               }
           }
           assemblyRegion = null;

           if (monitor.tooManyEnginesForLocalAssembly && monitor.localAssemblyRunnableLock.tryLock()) {
               if (monitor.tooManyEnginesForLocalAssembly) {
                   monitor.tooManyEnginesForLocalAssembly = false;
                   monitor.localAssemblyRunnableLock.unlock();
                   break;
               }
               monitor.localAssemblyRunnableLock.unlock();
           }
       }

       monitor.localAssemblyRunnableLock.lock();
       monitor.numRemainedThreadsForLocalAssembly.decrementAndGet();
       monitor.localAssemblyRunnableLock.unlock();
    }

    private AssemblyResultSet doWork (final AssemblyRegion originalAssemblyRegion, final FeatureContext featureContext) {
        final List<VariantContext> givenAlleles = MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ?
                featureContext.getValues(MTAC.alleles).stream().filter(vc -> MTAC.genotypeFilteredAlleles || vc.isNotFiltered()).collect(Collectors.toList()) :
                Collections.emptyList();

        final AssemblyRegion assemblyActiveRegion = AssemblyBasedCallerUtils.assemblyRegionWithWellMappedReads(originalAssemblyRegion, READ_QUALITY_FILTER_THRESHOLD, header);
        final AssemblyResultSet untrimmedAssemblyResult = AssemblyBasedCallerUtils.assembleReads(assemblyActiveRegion, givenAlleles, MTAC, header, samplesList, logger, referenceReader, assemblyEngine, aligner);
        final SortedSet<VariantContext> allVariantEvents = untrimmedAssemblyResult.getVariationEvents(MTAC.maxMnpDistance);
        final AssemblyRegionTrimmer.Result trimmingResult = trimmer.trim(originalAssemblyRegion, allVariantEvents);

        if (!trimmingResult.isVariationPresent()) {
            return null;
        }

        final AssemblyResultSet assemblyResult = trimmingResult.needsTrimming() ?
                untrimmedAssemblyResult.trimTo(trimmingResult.getCallableRegion()) : untrimmedAssemblyResult;

        if (!assemblyResult.isVariationPresent()) {
            return null;
        }

        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        final List<GATKRead> readStubs = regionForGenotyping.getReads().stream().filter(r -> r.getLength()<MINIMUM_READ_LENGTH_AFTER_TRIMMING).collect(Collectors.toList());
        regionForGenotyping.removeAll(readStubs);

        assemblyResult.reads = AssemblyBasedCallerUtils.splitReadsBySample(samplesList, header, regionForGenotyping.getReads());
        assemblyResult.extendedSpan = originalAssemblyRegion.getExtendedSpan();
        assemblyResult.tid = header.getSequenceIndex(originalAssemblyRegion.getContig());
        assemblyResult.regionStart = originalAssemblyRegion.getStart();
        assemblyResult.regionEnd = originalAssemblyRegion.getEnd();
        assemblyResult.givenAlleles = givenAlleles;

        return assemblyResult;
    }
}
