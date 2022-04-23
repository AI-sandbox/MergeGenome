## DataFix

When merging DNA sequences from different sources, it is important to remove any SNP that, by itself or combined, allows identifying the source of the sequences. The SNP distribution of zeros and ones should be similar for all sources with the same species and ancestries. However, no tool encompasses filtering such – likely erroneous – SNPs.

To identify and remove the features that, due to a sequencing error or a bad imputation, contribute most to identifying the source of the samples, the [DataFix](https://github.com/AI-sandbox/Datafix) software is of great utility. It uses ML models to detect and/or fix mismatching features between the reference and the query. Such as when imputing with Beagle, the former is assumed to contain high-quality data, while the latter might contain some - partially or completely - noisy features.

