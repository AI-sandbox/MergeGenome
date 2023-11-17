## SNP Error Correction

When erroneous SNPs are detected by the [binary classifiers](https://github.com/AI-sandbox/MergeGenome/blob/main/readmes/README_6_ml_source_classifiers_for_snp_filtering.md) trained for SNP filtering, the easiest approach is to remove them. However, this can lead to a substantial loss of valuable information. In many applications, it can be preferable to instead employ techniques that correct the erroneous SNPs. During the final stage of the MergeGenome's homogenization process, a multi-output machine learning classifier is trained to predict better SNP values to replace the filtered erroneous SNPs within the query dataset. According to our experiments, K-Nearest Neighbors (KNN) is the most suitable choice to correct erroneous SNPs caused by poor Beagle imputation and homogenization errors.

A newer version of the code can be found [here](https://github.com/AI-sandbox/Datafix).
