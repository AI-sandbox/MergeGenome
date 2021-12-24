#############################################
### Process genetic map for beagle
#############################################

import pandas as pd

# process genetic map for beagle

for chm in range(1,39):
    genetic_in = "/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/genetic_map/chr{}_average_canFam3.1.txt".format(chm)
    data = pd.read_csv(genetic_in, delimiter="\t", header="infer", comment="#",dtype=float)

    data = data.drop(data.columns[2],axis=1)

    genetic_out = "/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing/genetic_map/formatted/chr{0}_beagle_gmap.txt".format(chm)

    with open(genetic_out,"w") as f:
        for itr in data.iterrows():
            i=itr[1]
            #print(i[0],i[1],i[2])
            f.write("{}\t{}\t{}\t{}\n".format("chr"+str(int(i[0])),".",i[2],int(i[1])))