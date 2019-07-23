#!/usr/bin/python

import os
import re
import sys
import shutil

re_import = re.compile (r"^import")
re_tail = re.compile (r"^}")

if (len(sys.argv) < 3):
  print ("Usage: python update.py <update.dest.dir> <update.package>")
  sys.exit (-1)

gatk = os.path.abspath (sys.argv[1])
up_pack  = os.path.abspath (sys.argv[2])

#---------------------------------------------------------#
#                     Added Classes                       #
#---------------------------------------------------------#

if (os.path.exists(gatk + "/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Acc")):
  shutil.rmtree (gatk + "/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Acc")
shutil.copytree (up_pack+"/added_classes/Mutect2Acc", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Acc")

#---------------------------------------------------------#
#                    Modified Classes                     #
#---------------------------------------------------------#

###### Modify AssemblyRegionWalker

changed_lines = set ((
      "private List<MultiIntervalLocalReadShard> readShards;",
      "final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? userIntervals : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());",
      "readShards = makeReadShards(intervals);",
      "public final void traverse() {"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java", "r") as fp:
  for line in fp:
    info = line.strip ()
    if (info == "private List<MultiIntervalLocalReadShard> readShards;"):
      new_src.write ("    public List<MultiIntervalLocalReadShard> readShards;\n")
    elif (info == "final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? userIntervals : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());"):
      new_src.write ("        //" + info + "\n")
    elif (info == "readShards = makeReadShards(intervals);"):
      new_src.write ("        //" + info + "\n")
    elif (info == "public final void traverse() {"):
      new_src.write ("    public void traverse() {\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java")

###### Modify GATKTool

changed_lines = set ((
      "ReferenceDataSource reference;",
      "ReadsDataSource reads;",
      "FeatureManager features;",
      "List<SimpleInterval> userIntervals;"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/engine/GATKTool.java") as fp:
  for line in fp:
    info = line.strip ()
    if (info == "ReferenceDataSource reference;"):
      new_src.write ("    protected ReferenceDataSource reference;\n")
    elif (info == "ReadsDataSource reads;"):
      new_src.write ("    public ReadsDataSource reads;\n")
    elif (info == "FeatureManager features;"):
      new_src.write ("    protected FeatureManager features;\n")
    elif (info == "List<SimpleInterval> userIntervals;"):
      new_src.write ("    public List<SimpleInterval> userIntervals;\n");
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/engine/GATKTool.java")

###### Modify MultiIntervalLocalReadShard

import_added = False
func1_added = False
changed_lines = set ((
      "private final ReadsDataSource readsSource;"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/engine/MultiIntervalLocalReadShard.java") as fp:
  for line in fp:
    info = line.strip ()
    if (info == "private final ReadsDataSource readsSource;"):
      new_src.write ("    private ReadsDataSource readsSource;\n")
    elif ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import htsjdk.samtools.SAMSequenceDictionary;\n")
      new_src.write (line)
    elif (info == "public MultiIntervalLocalReadShard(final List<SimpleInterval> intervals, final int intervalPadding, final ReadsDataSource readsSource) {"):
      new_src.write ("    public MultiIntervalLocalReadShard(final List<SimpleInterval> intervals, final int intervalPadding, final SAMSequenceDictionary sequenceDictionary) {\n")
      new_src.write ("        Utils.nonNull(intervals);\n")
      new_src.write ('        Utils.validateArg(intervalPadding >= 0, "intervalPadding must be >= 0");\n')
      new_src.write ("\n")
      new_src.write ("        // Feed intervals through IntervalUtils.getIntervalsWithFlanks() to ensure they get sorted using\n")
      new_src.write ("        // the same comparator as the paddedIntervals below.\n")
      new_src.write ("        this.intervals = Collections.unmodifiableList(IntervalUtils.getIntervalsWithFlanks(intervals, 0, sequenceDictionary));\n")
      new_src.write ("\n")
      new_src.write ("        // This will both pad each interval and merge any intervals that are overlapping or adjacent after padding,\n")
      new_src.write ("        // in addition to sorting the intervals\n")
      new_src.write ("        this.paddedIntervals = Collections.unmodifiableList(IntervalUtils.getIntervalsWithFlanks(intervals, intervalPadding, sequenceDictionary));\n")
      new_src.write ("    }")
      new_src.write ("\n")
      new_src.write (line)
    elif ((not func1_added) and (info == "@Override")):
      func1_added = True
      new_src.write ("    public void setReadsSource (final ReadsDataSource readsSource) {\n")
      new_src.write ("        this.readsSource = readsSource;\n")
      new_src.write ("    }\n")
      new_src.write ("\n")
      new_src.write (line)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/engine/MultiIntervalLocalReadShard.java")

###### Modify AssemblyResultSet

import_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyResultSet.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.utils.read.GATKRead;\n")
      new_src.write ("import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;\n")
      new_src.write (line)
    elif (info == "public final class AssemblyResultSet {"):
      new_src.write (line)
      new_src.write ("\n")
      new_src.write ("    public Map<String,List<GATKRead>> reads;\n")
      new_src.write ("    public ReadLikelihoods<Haplotype> readLikelihoods;\n")
      new_src.write ("    public List<VariantContext> givenAlleles;\n")
      new_src.write ("    public SimpleInterval extendedSpan;\n")
      new_src.write ("    public int tid;\n")
      new_src.write ("    public int regionStart;\n")
      new_src.write ("    public int regionEnd;\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyResultSet.java")

###### Modify PairHMMLikelihoodCalculationEngine

import_added = False
func1_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMLikelihoodCalculationEngine.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.utils.pairhmm.LoglessPairHMM;\n")
      new_src.write ("import org.broadinstitute.hellbender.utils.pairhmm.VectorLoglessPairHMM;\n")
      new_src.write (line)
    elif (info == "final byte baseQualityScoreThreshold) {"):
      new_src.write (line)
      new_src.write ("        this (\n")
      new_src.write ("                constantGCP,\n")
      new_src.write ("                arguments,\n")
      new_src.write ("                hmmType,\n")
      new_src.write ("                log10globalReadMismappingRate,\n")
      new_src.write ("                pcrErrorModel,\n")
      new_src.write ("                baseQualityScoreThreshold,\n")
      new_src.write ("                false\n")
      new_src.write ("        );")
      new_src.write ("    }")
      new_src.write ("\n")
      new_src.write ("    public PairHMMLikelihoodCalculationEngine(final byte constantGCP,\n")
      new_src.write ("                                              final PairHMMNativeArguments arguments,\n")
      new_src.write ("                                              final PairHMM.Implementation hmmType,\n")
      new_src.write ("                                              final double log10globalReadMismappingRate,\n")
      new_src.write ("                                              final PCRErrorModel pcrErrorModel,\n")
      new_src.write ("                                              final byte baseQualityScoreThreshold,\n")
      new_src.write ("                                              final boolean onlyCPUPairHMMEngine) {\n")
    elif (info == "this.pairHMM = hmmType.makeNewHMM(arguments);"):
      new_src.write ("        if (onlyCPUPairHMMEngine) {\n")
      new_src.write ("            this.pairHMM = new LoglessPairHMM();\n")
      new_src.write ("        } else {\n")
      new_src.write ("            this.pairHMM = hmmType.makeNewHMM(arguments);\n")
      new_src.write ("        }\n")
    elif ((not func1_added) and (re_tail.search(line))):
      func1_added = True
      new_src.write ("\n")
      new_src.write ("    @Override\n")
      new_src.write ("    public VectorLoglessPairHMM.Implementation getRunningEngineType() {\n")
      new_src.write ("        return pairHMM.runningEngineType;\n")
      new_src.write ("    }\n")
      new_src.write (line)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMLikelihoodCalculationEngine.java")

###### Modify PairHMMNativeArgumentCollection

import_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMNativeArgumentCollection.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc.Mutect2Acc;\n")
      new_src.write (line)
    elif (info == '@Argument(fullName = "native-pair-hmm-threads", doc="How many threads should a native pairHMM implementation use", optional = true)'):
      new_src.write ('    //@Argument(fullName = "native-pair-hmm-threads", doc="How many threads should a native pairHMM implementation use", optional = true)\n')
    elif (info == "private int pairHmmNativeThreads = 4;"):
      new_src.write ("    //private int pairHmmNativeThreads = 4;\n")
    elif (info == "args.maxNumberOfThreads = pairHmmNativeThreads;"):
      new_src.write ("        if (Mutect2Acc.pairHmmNativeThreads <= 0) {\n")
      new_src.write ("            Mutect2Acc.pairHmmNativeThreads = Runtime.getRuntime().availableProcessors() / 2;\n")
      new_src.write ("        }\n")
      new_src.write ("        args.maxNumberOfThreads = Mutect2Acc.pairHmmNativeThreads;\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMNativeArgumentCollection.java")

###### Modify RandomLikelihoodCalculationEngine

import_added = False
func1_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/RandomLikelihoodCalculationEngine.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;\n")
      new_src.write ("import org.broadinstitute.hellbender.utils.pairhmm.VectorLoglessPairHMM;\n")
      new_src.write (line)
    elif ((not func1_added) and (re_tail.search(line))):
      func1_added = True
      new_src.write ("\n")
      new_src.write ("    @Override\n")
      new_src.write ("    public VectorLoglessPairHMM.Implementation getRunningEngineType() {\n")
      new_src.write ("        return null;\n")
      new_src.write ("    }\n")
      new_src.write (line)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/RandomLikelihoodCalculationEngine.java")

###### Modify ReadLikelihoodCalculationEngine

import_added = False
func1_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReadLikelihoodCalculationEngine.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.utils.pairhmm.VectorLoglessPairHMM;\n")
      new_src.write (line)
    elif ((not func1_added) and (re_tail.search(line))):
      func1_added = True
      new_src.write ("\n")
      new_src.write ("    public VectorLoglessPairHMM.Implementation getRunningEngineType();\n")
      new_src.write (line)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReadLikelihoodCalculationEngine.java")

###### Modify M2ArgumentCollection

changed_lines = set ((
      "protected String tumorSample = null;",
      "protected String normalSample = null;"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/M2ArgumentCollection.java") as fp:
  for line in fp:
    info = line.strip ()
    if (info == "protected String tumorSample = null;"):
      new_src.write ("    public String tumorSample = null;\n")
    elif (info == "protected String normalSample = null;"):
      new_src.write ("    public String normalSample = null;\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/M2ArgumentCollection.java")

###### Modify Mutect2Engine

changed_lines = set ((
      'private static final String MUTECT_VERSION = "2.1";',
      "private static List<Byte> altQuals(final ReadPileup pileup, final byte refBase, final int pcrErrorQual) {",
      "private static double lnLikelihoodRatio(final int refCount, final List<Byte> altQuals) {"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Engine.java") as fp:
  for line in fp:
    info = line.strip ()
    if (info == 'private static final String MUTECT_VERSION = "2.1";'):
      new_src.write ('    public static final String MUTECT_VERSION = "2.1";\n')
    elif (info == "private static List<Byte> altQuals(final ReadPileup pileup, final byte refBase, final int pcrErrorQual) {"):
      new_src.write ("    public static List<Byte> altQuals(final ReadPileup pileup, final byte refBase, final int pcrErrorQual) {\n")
    elif (info == "private static double lnLikelihoodRatio(final int refCount, final List<Byte> altQuals) {"):
      new_src.write ("    public static double lnLikelihoodRatio(final int refCount, final List<Byte> altQuals) {")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Engine.java")

###### Modify GenotypeUtils

changed_lines = set ((
      "final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2);",
      "final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a1);"))

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/utils/GenotypeUtils.java") as fp:
  for line in fp:
    info = line.strip ()
    if (info == "final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2);"):
      new_src.write ("                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(a2);\n")
    elif (info == "final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a1);"):
      new_src.write ("                    final int[] idxVector = vc.getGLIndicesOfAlternateAllele(a1);\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/utils/GenotypeUtils.java")

###### Modify PairHMM

import_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not import_added) and (re_import.search(info))):
      import_added = True
      new_src.write ("import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Acc.Mutect2Acc;\n")
      new_src.write (line)
    elif (info == "public abstract class PairHMM implements Closeable{"):
      new_src.write (line)
      new_src.write ("\n")
      new_src.write ("    public VectorLoglessPairHMM.Implementation runningEngineType;\n")
    elif (info == "FASTEST_AVAILABLE(args -> {"):
      new_src.write (line)
      new_src.write ("\n")
      new_src.write ("            if (Mutect2Acc.useFPGA) {\n")
      new_src.write ("                try {\n")
      new_src.write ("                    final VectorLoglessPairHMM hmm = new VectorLoglessPairHMM(VectorLoglessPairHMM.Implementation.FPGA, args);\n")
      new_src.write ('                    logger.info("Using the FPGA-accelerated native PairHMM implementation");\n')
      new_src.write ("                    return hmm;\n")
      new_src.write ("                } catch (UserException.HardwareFeatureException e) {\n")
      new_src.write ('                    logger.info("FPGA-accelerated native PairHMM implementation is not supported");\n')
      new_src.write ("                }\n")
      new_src.write ("            }\n")
      new_src.write ("\n")
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java")

###### Modify VectorLoglessPairHMM

n_skip_lines = 0

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/VectorLoglessPairHMM.java") as fp:
  for line in fp:
    if (n_skip_lines > 0):
      n_skip_lines -= 1
      continue

    info = line.strip ()
    if (info == "public final class VectorLoglessPairHMM extends LoglessPairHMM {"):
      new_src.write (line)
      new_src.write ("\n")
      new_src.write ("    private int NUM_READS_IN_BATCH = 16384;\n")
    elif (info == 'throw new UserException.HardwareFeatureException("Machine does not support AVX PairHMM.");'):
      new_src.write (line)
      new_src.write ("                }\n")
      new_src.write ("                runningEngineType = Implementation.AVX;\n")
      n_skip_lines = 1
    elif (info == 'throw new UserException.HardwareFeatureException("Machine does not support OpenMP AVX PairHMM.");'):
      new_src.write (line)
      new_src.write ("                }\n")
      new_src.write ("                runningEngineType = Implementation.OMP;\n")
      n_skip_lines = 1
    elif (info == 'throw new UserException.HardwareFeatureException("Machine does not support FPGA PairHMM.");'):
      new_src.write (line)
      new_src.write ("                }\n")
      new_src.write ("                runningEngineType = Implementation.FPGA;\n")
      n_skip_lines = 1
    elif (info == "ReadDataHolder[] readDataArray = new ReadDataHolder[readListSize];"):
      new_src.write ("        //ReadDataHolder[] readDataArray = new ReadDataHolder[readListSize];\n")
      new_src.write ("        //int idx = 0;\n")
      new_src.write ("        //for (GATKRead read : processedReads) {\n")
      new_src.write ("        //    readDataArray[idx] = new ReadDataHolder();\n")
      new_src.write ("        //    readDataArray[idx].readBases = read.getBases();\n")
      new_src.write ("        //    readDataArray[idx].readQuals = read.getBaseQualities();\n")
      new_src.write ("        //    readDataArray[idx].insertionGOP = ReadUtils.getBaseInsertionQualities(read);\n")
      new_src.write ("        //    readDataArray[idx].deletionGOP = ReadUtils.getBaseDeletionQualities(read);\n")
      new_src.write ("        //    readDataArray[idx].overallGCP = gcp.get(read);\n")
      new_src.write ("        //    ++idx;\n")
      new_src.write ("        //}\n")
      n_skip_lines = 10
    elif (info == "pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, mLogLikelihoodArray);"):
      new_src.write ("        //pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, mLogLikelihoodArray);\n")
      new_src.write ("\n")
      code = '''
        if (readListSize>=32768 && runningEngineType==Implementation.FPGA) {
            int readIdx = 0;
            int batchIdx = 0;
            int numBatch = (readListSize + NUM_READS_IN_BATCH - 1) / NUM_READS_IN_BATCH;
            int begIndexForResult = 0;
            int numOneBatchLikelihood = NUM_READS_IN_BATCH * numHaplotypes;
            double[] likelihoodArray;
            ReadDataHolder[] readDataArray;

            int numBatchRead;
            int numBatchLikelihood;
            if (batchIdx == numBatch-1) {
                numBatchRead = readListSize;
                numBatchLikelihood = readListSize * numHaplotypes;
            } else {
                numBatchRead = NUM_READS_IN_BATCH;
                numBatchLikelihood = numOneBatchLikelihood;
            }
            readDataArray = new ReadDataHolder[numBatchRead];
            likelihoodArray = new double[numBatchLikelihood];
            for (int i=0; i<numBatchRead; ++i) {
                readDataArray[i] = new ReadDataHolder();
            }

            for (GATKRead read : processedReads) {
                readDataArray[readIdx].readBases = read.getBases();
                readDataArray[readIdx].readQuals = read.getBaseQualities();
                readDataArray[readIdx].insertionGOP = ReadUtils.getBaseInsertionQualities(read);
                readDataArray[readIdx].deletionGOP = ReadUtils.getBaseDeletionQualities(read);
                readDataArray[readIdx].overallGCP = gcp.get(read);
                ++readIdx;
                if (readIdx == NUM_READS_IN_BATCH) {
                    readIdx = 0;
                    pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, likelihoodArray);
                    System.arraycopy(likelihoodArray, 0, mLogLikelihoodArray, begIndexForResult, numBatchLikelihood);
                    begIndexForResult += numBatchLikelihood;
                    ++batchIdx;

                    if (batchIdx == numBatch-1) {
                        numBatchRead = readListSize - batchIdx*NUM_READS_IN_BATCH;
                        if (numBatchRead <= 0) {
                            break;
                        }
                        numBatchLikelihood = numBatchRead * numHaplotypes;
                    } else {
                        numBatchRead = NUM_READS_IN_BATCH;
                        numBatchLikelihood = numOneBatchLikelihood;
                    }
                    readDataArray = new ReadDataHolder[numBatchRead];
                    for (int i=0; i<numBatchRead; ++i) {
                        readDataArray[i] = new ReadDataHolder();
                    }
                    likelihoodArray = new double[numBatchLikelihood];
                }
            }
            if (readIdx != NUM_READS_IN_BATCH) {
                pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, likelihoodArray);
                System.arraycopy(likelihoodArray, 0, mLogLikelihoodArray, begIndexForResult, numBatchLikelihood);
            }
        } else {
            ReadDataHolder[] readDataArray = new ReadDataHolder[readListSize];
            int idx = 0;
            for (GATKRead read : processedReads) {
                readDataArray[idx] = new ReadDataHolder();
                readDataArray[idx].readBases = read.getBases();
                readDataArray[idx].readQuals = read.getBaseQualities();
                readDataArray[idx].insertionGOP = ReadUtils.getBaseInsertionQualities(read);
                readDataArray[idx].deletionGOP = ReadUtils.getBaseDeletionQualities(read);
                readDataArray[idx].overallGCP = gcp.get(read);
                ++idx;
            }

            //for(reads)
            //   for(haplotypes)
            //       compute_full_prob()
            //long beg = System.nanoTime();
            pairHmm.computeLikelihoods(readDataArray, mHaplotypeDataArray, mLogLikelihoodArray);
            //System.out.printf ("\\t\\t\\t\\tvector computeReadLikelihoods cost: %fs\\n", (double)(System.nanoTime()-beg) / 1000000000.0);

        }
      '''
      new_src.write (code)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/VectorLoglessPairHMM.java")

###### Modify GVCFWriter

func1_added = False

new_src = open ("tmp.java", "w")
with open (gatk+"/src/main/java/org/broadinstitute/hellbender/utils/variant/writers/GVCFWriter.java") as fp:
  for line in fp:
    info = line.strip ()
    if ((not func1_added) and (info == "@Override")):
      func1_added = True
      new_src.write ("    @Override\n")
      new_src.write ("    public void dump(VariantContext vc) {\n")
      new_src.write ("\n")
      new_src.write ("    }\n")
      new_src.write ("\n")
      new_src.write (line)
    else:
      new_src.write (line)
new_src.close ()

shutil.copy ("tmp.java", gatk+"/src/main/java/org/broadinstitute/hellbender/utils/variant/writers/GVCFWriter.java")
