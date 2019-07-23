package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticGenotypingEngine;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class GenotypingRunnable implements Runnable {
    private final MonitorRunnable monitor;
    private final CommandLineProgram instance;
    private final M2ArgumentCollection MTAC;
    private final ReferenceInputArgumentCollection referenceArguments;
    private final int cloudPrefetchBuffer;
    private final int cloudIndexPrefetchBuffer;
    private final SampleList samplesList;
    private final String tumorSample;
    private final String normalSample;
    private final SAMFileHeader header;

    private ReferenceDataSource reference;
    private FeatureManager features;
    private SmithWatermanAligner aligner;
    private VariantAnnotatorEngine annotatorEngine;
    private SomaticGenotypingEngine genotypingEngine;
    private VCFEncoder vcfEncoder;

    GenotypingRunnable (
            final MonitorRunnable monitor,
            final CommandLineProgram instance,
            final M2ArgumentCollection MTAC,
            final ReferenceInputArgumentCollection referenceArguments,
            final int cloudPrefetchBuffer,
            final int cloudIndexPrefetchBuffer,
            final SampleList samplesList,
            final String tumorSample,
            final String normalSample,
            final Collection<Annotation> annotations,
            final SAMFileHeader header,
            final VCFHeader vcfHeader,
            final boolean allowMissingFieldsInHeader
    ) {
        this.monitor = monitor;
        this.instance = instance;
        this.MTAC = MTAC;
        this.referenceArguments = referenceArguments;
        this.cloudPrefetchBuffer = cloudPrefetchBuffer;
        this.cloudIndexPrefetchBuffer = cloudIndexPrefetchBuffer;
        this.samplesList = samplesList;
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
        this.header = header;

        reference = Mutect2Acc.LocalInitializeReference(referenceArguments);
        features = Mutect2Acc.LocalInitializeFeatures(instance, referenceArguments, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        aligner = SmithWatermanAligner.getAligner(MTAC.smithWatermanImplementation);
        annotatorEngine = new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false);
        genotypingEngine = new SomaticGenotypingEngine(samplesList, MTAC, tumorSample, normalSample);
        genotypingEngine.setAnnotationEngine(annotatorEngine);
        vcfEncoder = new VCFEncoder(vcfHeader, allowMissingFieldsInHeader, false);
    }

    public void run () {
        monitor.numInitializedThreadsForGenotyping.incrementAndGet();

        monitor.genotypingRunnableLock.lock();
        monitor.numRemainedThreadsForGenotyping.incrementAndGet();
        monitor.genotypingRunnableLock.unlock();

        monitor.runnableIsStarted = true;

        while (true) {
            AssemblyResultSet assemblyResult = monitor.pairHMMResultQueue.poll();
            if (assemblyResult == null) {
                if (monitor.vectorLogessPairHMMWorkDone && monitor.numRemainedThreadsForLoglessPairHMM.get()<=0) {
                    break;
                }
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
                continue;
            }

            List<VariantContext> variants = doWork (
                    assemblyResult,
                    new ReferenceContext(reference, assemblyResult.extendedSpan),
                    new FeatureContext(features, assemblyResult.extendedSpan)
            );

            monitor.numActiveRegionAddressedByGenotyping.incrementAndGet();

            if (variants == null) {
                continue;
            }

            ActiveRegionVariantSet set = new ActiveRegionVariantSet(
                    variants,
                    assemblyResult.tid,
                    assemblyResult.regionStart,
                    assemblyResult.regionEnd
            );

            while (!monitor.variantSetQueue.offer(set)) {
                try {
                    Thread.sleep(1);
                } catch (InterruptedException e) {

                }
            }
            assemblyResult = null;

            if (monitor.tooManyEnginesForGenotyping && monitor.genotypingRunnableLock.tryLock()) {
                if (monitor.tooManyEnginesForGenotyping) {
                    monitor.tooManyEnginesForGenotyping = false;
                    monitor.genotypingRunnableLock.unlock();
                    break;
                }
                monitor.genotypingRunnableLock.unlock();
            }
        }

        monitor.genotypingRunnableLock.lock();
        monitor.numRemainedThreadsForGenotyping.decrementAndGet();
        monitor.genotypingRunnableLock.unlock();
    }

    private List<VariantContext> doWork (
            final AssemblyResultSet assemblyResult,
            final ReferenceContext referenceContext,
            final FeatureContext featureContext
    ) {
        ReadLikelihoods<Haplotype> readLikelihoods = assemblyResult.readLikelihoods;
        final AssemblyRegion regionForGenotyping = assemblyResult.getRegionForGenotyping();
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(
                readLikelihoods,
                assemblyResult.getReferenceHaplotype(),
                assemblyResult.getPaddedReferenceLoc(),
                aligner
        );
        readLikelihoods.changeReads(readRealignments);

        final HaplotypeCallerGenotypingEngine.CalledHaplotypes calledHaplotypes = genotypingEngine.callMutations(
                readLikelihoods,
                assemblyResult,
                referenceContext,
                regionForGenotyping.getSpan(),
                featureContext,
                assemblyResult.givenAlleles,
                header
        );

        List<VariantContext> variants = calledHaplotypes.getCalls();

        if (variants.size() <= 0) {
            return null;
        }

        for (final VariantContext var : variants) {
            var.dumpString = vcfEncoder.encode(var);
        }

        return variants;
    }
}
