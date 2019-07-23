package org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

class ActiveRegionVariantSet implements Comparable<ActiveRegionVariantSet> {
    final public List<VariantContext> vaiants;
    final int tid;
    final int regionStart;
    final int regionEnd;

    ActiveRegionVariantSet (
            final List<VariantContext> vaiants,
            final int tid,
            final int regionStart,
            final int regionEnd) {
        this.vaiants = vaiants;
        this.tid = tid;
        this.regionStart = regionStart;
        this.regionEnd = regionEnd;
    }

    @Override
    public int compareTo (ActiveRegionVariantSet s) {
        if (this.tid > s.tid)
            return 1;
        if (this.tid < s.tid)
            return -1;

        return (this.regionStart - s.regionStart);
    }
}
