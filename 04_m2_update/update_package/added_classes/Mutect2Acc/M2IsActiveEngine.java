package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingGivenAllelesUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingOutputMode;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.util.List;

public class M2IsActiveEngine implements AssemblyRegionEvaluator {
    final M2ArgumentCollection MTAC;
    final SAMFileHeader header;
    final Logger logger;
    final String tumorSample;
    final String normalSample;

    public static final double MAX_ALT_FRACTION_IN_NORMAL = 0.3;
    public static final int MAX_NORMAL_QUAL_SUM = 100;

    public M2IsActiveEngine (final M2ArgumentCollection MTAC,
                             final SAMFileHeader header,
                             final Logger logger,
                             final String tumorSample,
                             final String normalSample) {
        this.MTAC = MTAC;
        this.header = header;
        this.logger = logger;
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
    }

    @Override
    public ActivityProfileState isActive (final AlignmentContext context,
                                          final ReferenceContext ref,
                                          final FeatureContext featureContext) {
        final SimpleInterval refInterval = ref.getInterval();
        if (MTAC.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES) {
            final VariantContext vcFromAllelesRod = GenotypingGivenAllelesUtils.composeGivenAllelesVariantContextFromRod(featureContext, refInterval, false, MTAC.genotypeFilteredAlleles, logger, MTAC.alleles);
            if (vcFromAllelesRod != null) {
                return new ActivityProfileState(refInterval, 1.0);
            }
        }

        final byte refBase = ref.getBase();
        if (context==null || context.getBasePileup().isEmpty()) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        final ReadPileup pileup = context.getBasePileup();
        final ReadPileup tumorPileup = pileup.getPileupForSample(tumorSample, header);
        final List<Byte> tumorAltQuals = Mutect2Engine.altQuals(tumorPileup, refBase);
        final double tumorLog10Odds = MathUtils.logToLog10(Mutect2Engine.lnLikelihoodRatio(tumorPileup.size()-tumorAltQuals.size(), tumorAltQuals));

        if (tumorLog10Odds < MTAC.initialTumorLod) {
            return new ActivityProfileState(refInterval, 0.0);
        } else if (normalSample!=null && !MTAC.genotypeGermlineSites) {
            final ReadPileup normalPileup = pileup.getPileupForSample(normalSample, header);
            final List<Byte> normalAltQuals = Mutect2Engine.altQuals(normalPileup, refBase);
            final int normalAltCount = normalAltQuals.size();
            final double normalQualSum = normalAltQuals.stream().mapToDouble(Byte::doubleValue).sum();
            if (normalAltCount>normalPileup.size()*MAX_ALT_FRACTION_IN_NORMAL && normalQualSum>MAX_NORMAL_QUAL_SUM) {
                return new ActivityProfileState(refInterval, 0.0);
            }
        } else if (!MTAC.genotypeGermlineSites) {
            final List<VariantContext> germline = featureContext.getValues(MTAC.germlineResource, refInterval);
            if (!germline.isEmpty()) {
                final List<Double> germlineAlleleFrequnecies = germline.get(0).getAttributeAsDoubleList(VCFConstants.ALLELE_COUNT_KEY, 0.0);
                if (!germlineAlleleFrequnecies.isEmpty()
                        && germlineAlleleFrequnecies.get(0) > MTAC.maxPopulationAlleleFrequency) {
                    return new ActivityProfileState(refInterval, 0.0);
                }
            }
        }

        if (!MTAC.genotypePonSites && !featureContext.getValues(MTAC.pon, new SimpleInterval(context.getContig(), (int)context.getPosition(), (int)context.getPosition())).isEmpty()) {
            return new ActivityProfileState(refInterval, 0.0);
        }

        return new ActivityProfileState(refInterval, 1.0, ActivityProfileState.Type.NONE, null);
    }
}
