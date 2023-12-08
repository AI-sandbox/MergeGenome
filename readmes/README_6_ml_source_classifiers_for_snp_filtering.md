## Machine Learning Source Classifiers for SNP filtering

When combining DNA sequences from two different sources, it's critical to eliminate any SNP that, by itself or combined, reveals the source of the samples. As far as the DNA sequences are from similar species or ancestries, the distribution of zeros and ones should be similar between SNPs coming from different origins. In case there is a substantial difference in the SNP distributions, a sequencing error most certainly occurred. Alternatively, in the situation when missing allele data was inferred through an imputation model, a change in the SNP means between sources could suggest that the model used for inferring a portion of the allele data was insufficiently accurate. Such SNPs might irreversibly contaminate any analysis constructed on top of it, so they should be carefully dealt with.

A newer version of this functionality can be found [here](https://github.com/AI-sandbox/Datafix).
