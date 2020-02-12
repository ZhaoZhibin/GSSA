# GSSA
Source codes for the paper "Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware"



This repository contains the implementation details of our paper: [TIM]
"[**Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware**](https://zhaozhibin.github.io/)" 
by [Zhibin Zhao](https://zhaozhibin.github.io/). 


## About
Vibration signal analysis has become one of the important methods for machinery fault diagnosis. Extraction of weak fault features from vibration signals with heavy background noise remains a challenging problem. In this paper, we first introduce the idea of algorithm-aware sparsity-assisted methods for fault feature enhancement, which extends model-aware sparsity-assisted fault diagnosis and allows a more flexible and convenient algorithm design. In the framework of algorithm-aware methods, we define the generalized structured shrinkage operators and construct the generalized structured shrinkage algorithm (GSSA) to overcome the disadvantages of $ l_1 $-norm regularization based fault feature enhancement methods. We then perform a series of simulation studies and two experimental cases to verify the effectiveness of the proposed method. Additionally, comparisons with model-aware methods, including basis pursuit denoising and windowed-group-lasso, and fast kurtogram further verify the advantages of GSSA for weak fault feature enhancement.

## Dependencies
- Matlab R2016b
- [[TQWT—toolbox]](http://eeweb.poly.edu/iselesni/TQWT/index.html) from I. W. Selesnick. 

## Pakages

This repository is organized as:
- [Generalized-Structured-Shrinkage-Operators](https://github.com/ZhaoZhibin/GSSA/tree/master/Generalized-Structured-Shrinkage-Operators) contains the main functions of the algorithm.
- [utils](https://github.com/ZhaoZhibin/GSSA/tree/master/utils) contains the extra functions of the test.
- [Performance_Data](https://github.com/ZhaoZhibin/GSSA/tree/master/Performance_Data) contains the performance comparison of different algorithms.
- [results](https://github.com/ZhaoZhibin/GSSA/tree/master/results) contains the results of the algorithm.
- [Neighborhood_Thresholding](https://github.com/ZhaoZhibin/GSSA/tree/master/Neighborhood_Thresholding) contains the TQWT denoising methods using neighborhood thresholding.
- [tqwt_matlab_toolbox](https://github.com/ZhaoZhibin/GSSA/tree/master/tqwt_matlab_toolbox) contains the TQWT toolbox copied from I. W. Selesnick.
In our implementation, **Matlab R2016b** is used to perform all the experiments.

Main functions:
- [Plot_Pure_Noise_Comparison.m] performs denoising when the features are submerged by Gaussian noise.
- [Plot_Noise_Plus_Harmonic_Interference_Comparison.m] performs denoising when the features are submerged by Gaussian noise and harmonic interference.
- [Plot_Penalty_Shrinkage.m] plots the penalties and their corresponding shrinkages.
- [Plot_performance_comparison.m] performs performance comparison of different algorithms.




## Implementation:
Flow the steps presented below:
-  Clone this repository.
```
git clone https://github.com/ZhaoZhibin/GSSA.git
open it with matlab
```
-  Test Simulation with Gaussian noise: Run `Plot_Pure_Noise_Comparison.m`. 
-  Test Simulation with Gaussian noise and harmonic interference: Run `Plot_Noise_Plus_Harmonic_Interference_Comparison.m`. 


## Citation
If you feel our GSSA is useful for your research, please consider citing our paper: 

```
@article{zhao2020sparsity,
  title={Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware},
  author={Zhao, Zhibin and Wang, Shibin and Xu, Weixin and Wu, Shuming and Wong David and Chen, Xuefeng},
  journal={IEEE Transactions on Instrumentation and Measurement},
  year={2020},
  publisher={IEEE}
}
```
## Contact
- zhibinzhao1993@gmail.com

