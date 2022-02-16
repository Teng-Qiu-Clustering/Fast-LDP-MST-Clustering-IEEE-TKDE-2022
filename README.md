# Fast LDP-MST: an efficient density-based clustering method for large-size datasets 

>**Teng Qiu, Yongjie Li, IEEE Transactions on Knowledge and Data Engineering, 2022, DOI: 10.1109/TKDE.2022.3150403.**

**Introduction**: Cluster analysis is widely studied and used in diverse
fields. Despite that
many clustering methods have been proposed, it is rare for a clustering
method to perform well simultaneously on the following characteristics:  
* effectiveness (in terms of accuracy);
* efficiency (in terms of speed);
* robustness (in terms of noise and parameter sensitivity);
* user-friendliness (in terms of the number of user-specified parameters,
interpretability, and reproducibility of the results, etc.). 

This work contributes to the field of clustering by providing a new solution (i.e., Fast LDP-MST) which largely improves the efficiency of LDP-MST, without sacrificing the merits of LDP-MST in effectiveness, robustness, and user-friendliness. 
Fast LDP-MST achieves a good balance among effectiveness, efficiency, robustness, and user-friendliness, and thus it could
have a certain degree of practical value in this big data era.

# Quick start
**Demo1.m** and **Demo2.m** were successfully tested on Matlab2018 (on a computer with a 3 GHz Intel i5 processor and 32 GB RAM). 

Note: **Besides downloading and running the codes, one can also directly run the code in my code ocean: https://codeocean.com/capsule/3297972/tree**

# Figures

![Fig.1](https://github.com/Teng-Qiu-Clustering/Fast-LDP-MST-Clustering/blob/main/Result/LDPMST-vs-FastLDPMST_on_TB_SF_CC_CG_Flower.png)
**Fig.1**: Comparison between LDP-MST and Fast LDP-MST on multiple datasets sampled from five different density functions (TB, SF, CC, CG, and Flower). The first and second rows compare the runtime (in seconds) and clustering accuracy
(ARI scores) of the two methods on a set of datasets with the numbers of samples varying from N=2^14 to N=2^24 (**> 16 million**). For instance, the test datasets with N=2^15 samples are displayed in the third row, where each point represents a sample; the corresponding clustering results of Fast LDP-MST are shown in the fourth row, where different colors on the points in each dataset indicate different clusters they are assigned to). The NMI scores and the visualization of the clustering results of LDP-MST are omitted here, as both methods obtain almost the same performance. ARI: adjusted Rand index (its value ranges
from -1 to 1; higher values
indicate higher clustering accuracy). 
 
 
 **** 
 
 ![Fig.2](https://github.com/Teng-Qiu-Clustering/Fast-LDP-MST-Clustering/blob/main/Result/LDPMST-vs-FastLDPMST_on_GSC.png)
**Fig.2**: Comparison between LDP-MST (2nd column) and Fast LDP-MST (3rd column)
on two challenging datasets (1st column). Each point represents a 2-dimensional sample. Different colors (in 2nd and 3rd columns) on points indicate different clustering labels assigned by the clustering methods. 
 
 
