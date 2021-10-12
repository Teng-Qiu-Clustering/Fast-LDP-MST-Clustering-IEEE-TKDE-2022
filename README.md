# FastLDPMST
 FastLDPMST: an efficient density-based clustering method for large-size datasets

# Quick start
**Demo1.m** and **Demo2.m** were successfully tested on Matlab2018

Note: 

1) for **Demo 1**, you can choose different datasets (**shown in the datasets section**) to test.

2) for **Demo 2**, the clustering results may be different in different runs, because the test datasets in Demo2 are randomly sampled from some density functions (i.e., the test datasets are different in different runs). Note that the proposed method is not sensitive to the initialization (i.e., for the same dataset, the proposed method will show the same result in different runs. This can be demonstrated by Demo 1). 

3) **You can also choose to directly run the code in my code ocean: https://codeocean.com/capsule/3297972/tree/v2**

# Figures

![Fig.1](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/LDPMST-vs-FastLDPMST_on_TB_SF_CC_CG_Flower.png)
**Fig.1**: Comparison of the efficiency of LDP-MST and FastLDPMST. The first and second rows compare the runtime (in seconds) and clustering
accuracy (ARI) of the two methods. The red arrows in the second row indicate the ARI scores of LDP-MST and FastLDPMST on the datasets with
N = 2^15, with the visualized clustering results of FastLDPMST being shown in the third row (the visualized clustering results of LDP-MST are
almost the same and thus omitted here). For each dataset in the third row, we use different colors to indicate different clusters that the data points
are assigned to. The results on NMI are omitted here, as both methods obtain almost the same performance.
  
 
****

![Fig.2](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/One_Dim_Uniform_V2.png)
**Fig.2**: Comparison between FastLDPMST and LDP-MST on 1-dimensional uniformly distributed datasets (with 2 clusters).

 
 **** 
 
 ![Fig.3](https://github.com/Teng-Qiu-Clustering/FastLDPMST/blob/main/LDPMST-vs-FastLDPMST_on_GSC.png)
**Fig.3**: Comparison between FastLDPMST and LDP-MST on two synthetic
datasets. The first row shows the test datasets, GaSpCi (left) and
GaSpCiNo (right). The second and third rows show the corresponding
clustering results obtained by FastLDPMST and LDP-MST, respectively.
 
