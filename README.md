# SJSONMF

The code in this toolbox implements ["Structured Joint Sparse Orthogonal Nonnegative Matrix Factorization for Fault Detection"](https://ieeexplore.ieee.org/document/10036023). More specifically, the abstract is given as follows. 

```As modern industrial processes become complicated, and some faults are difficult to be detected due to noises and nonlinearity of data, data-driven fault detection (FD) has been extensively used to detect abnormal events in functional units. To obtain better FD performance of nonnegative matrix factorization (NMF), this article first proposes an FD method using the structured joint sparse orthogonal NMF (SJSONMF). The core idea is to incorporate the graph regularization, sparsity, and orthogonality constraints into the classical NMF, which enjoys stronger discriminative ability, removes redundancy of different basis vectors, and improves fault interpretability. More importantly, an optimization algorithm based on the proximal alternating nonnegative least squares (PANLS) is developed, which can guarantee and speed up the convergence. Finally, the effectiveness of the proposed method is demonstrated by the experiments on the benchmark Tennessee Eastman Process (TEP) and two practical bearing datasets. Particularly, compared with the classical NMF, the T2 statistic has a gain of 33.13% for the fault IDV(16) on the TEP. The results show that the proposed model and algorithms are promising for FD.```




### Testing

Directly run demo_SJSONMF.m for reproduction.

### Citation
Please give credits to this paper if this code is useful and helpful for your research.

     @article{sun2023learning,
      title     = {Structured Joint Sparse Orthogonal Nonnegative Matrix Factorization for Fault Detection},
      author    = {Zhang, Xi and Xiu, Xianchao and Zhang, Chao},
      journal   = {IEEE Transactions on Instrumentation and Measurement},
      year      = {2023},
      volume    = {72},
      pages     = {1--15},
      publisher = {IEEE}
     }


### Contact 
Please feel free to contact xcxiu@shu.edu.cn if you have any questions.











