#!/usr/bin/python

out = open ("test.java", "w")

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
            //System.out.printf ("\t\t\t\tvector computeReadLikelihoods cost: %fs\\n", (double)(System.nanoTime()-beg) / 1000000000.0);

        }
'''

out.write (code)
out.close ()
