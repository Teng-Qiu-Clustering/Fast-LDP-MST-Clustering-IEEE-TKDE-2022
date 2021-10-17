# FastLDPMST: an efficient density-based clustering method for large-size datasets.
Cluster analysis is widely studied and used in diverse
fields. Despite that
many clustering methods have been proposed, it is rare for a clustering
method to perform well simultaneously on the following characteristics:  
* effectiveness (in terms of accuracy);
* efficiency (in terms of speed);
* robustness (in terms of noise and parameter sensitivity);
* user-friendliness (in terms of the number of user-specified parameters,
interpretability, and reproducibility of the results, etc.). 

This work contributes to the field of clustering by providing a new solution (i.e., FastLDPMST) which largely improves the efficiency of LDP-MST, without sacrificing the merits of LDP-MST in effectiveness, robustness, and user-friendliness. 
FastLDPMST achieves a good balance among effectiveness, efficiency, robustness, and user-friendliness, and thus it could
have a certain degree of practical value in this big data era.

# Quick start
**Demo1.m** and **Demo2.m** were successfully tested on Matlab2018 (on a computer with a 3 GHz Intel i5 processor and 32 GB RAM). 

Note: 

1) for **Demo 1**, one can choose different datasets (**shown in the datasets section**) to test.

2) for **Demo 2**, the clustering results may be different in different runs, because the test datasets in Demo2 are randomly sampled from some density functions (i.e., the test datasets are different in different runs). Note that the proposed method is not sensitive to the initialization (i.e., for the same dataset, the proposed method will show the same result in different runs. This can be demonstrated by Demo 1). 

3) **One can also choose to directly run the code in my code ocean: https://codeocean.com/capsule/3297972/tree/v2**

# Figures

![Fig.1](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/LDPMST-vs-FastLDPMST_on_TB_SF_CC_CG_Flower.png)
**Fig.1**: Comparison between LDP-MST and FastLDPMST on multiple datasets sampled from five different density functions (TB, SF, CC, CG, and Flower). The first and second rows compare the runtime (in seconds) and clustering accuracy
(ARI scores) of the two methods on a set of datasets with the numbers of samples varying from N=2^14 to N=2^24. For instance, the test datasets with N=2^15 samples are displayed in the third row, where each point represents a sample; the corresponding clustering results of FastLDPMST are shown in the fourth row, where different colors on the points in each dataset indicate different clusters they are assigned to). The NMI scores and the visualization of the clustering results of LDP-MST are omitted here, as both methods obtain almost the same performance. ARI: adjusted Rand index (its value ranges
from -1 to 1; higher values
indicate higher clustering accuracy). 
  
 
****

![Fig.2](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/One_Dim_Uniform_V2.png)
**Fig.2**: Comparison between FastLDPMST and LDP-MST on a set of 1-dimensional uniformly distributed datasets with the numbers of samples varying from N=2^14 to N=2^24. (a) shows how the ratios of the root nodes of the two methods vary with the number of samples, **which indicates that the number of root nodes could be as close as the number of samples in the worst case, which means that when testing such kinds of datasets, the runtime of LDPMST becomes O(N^2) in sharp contrast to ours (with O(N*logN))**. (b), (c) and (d) compare the runtimes of the two methods on three sub-steps. (e) compares the total runtimes (on all the eight steps) of the two methods. **Notably, as shown in (e), when N = 2^14, the total runtime of LDP-MST is 6495.2 seconds while that of FastLDPMST is 0.15 seconds, which means that the speedup factor is as high as 41512. Since the total runtime of LDP-MST increases much rapidly than that of FastLDPMST, the speedup factor would be much higher when N > 2^14.**

 
 **** 
 
 ![Fig.3](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/LDPMST-vs-FastLDPMST_on_GSC.png)
**Fig.3**: Comparison between LDP-MST (2nd column) and FastLDPMST (3rd column)
on two challenging datasets (1st column).
 
  **** 
**About the parameters:** Like LDP-MST, FastLDPMST contains two
parameters: K and eta. Parameter K is used as an early
termination condition for the number of nearest neighbors; Parameter
eta is used as a constraint for the cluster size. Note that in this
study, we do not regard cluster number as a parameter. For each dataset tested in all of the above figures, we fixed the two parameters k and eta to their emperical values (being log2(N) and 0.018*N, respectively) without careful tuning. In the following, we will show the influence of different values of k (with eta being fixed as its emperical value). 

 ![Fig.4](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/CompareFastLDPMST_and_DP_with_diff_parameters.png)
**Fig.4**: Comparison between FastLDPMST and DP (density-peak-based clustering) with different values of k (ranging from 7 to 100 with a step of 1). See examples in Fig. 5.

 ![Fig.5](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/Compare_visualization_of_clustering_Between_FastLDPMST_and_DP.png)
**Fig.5**: Comparison between FastLDPMST and DP (density-peak-based clustering) under the emperical value of k. For FastLDPMST, k = log2(N) (eta = 0.018*N); For DP, except dataset 3Circles (k = 30; there reported an error for DP while testing 3Circles with k = log2(N)), on all the other datasets, k = log2(N). 
