# Sparse Bayesian Multidimensional Item Response Theory

This repo contains the official R implementation of the PXL-EM algorithm to estimate sparse factor loadings from binary item response data, as illustrated in our paper *Sparse Bayesian Multidimensional Item Response Theory*: (https://arxiv.org/pdf/2310.17820.pdf). While we advocate for this Bayesian EM approach, we also publicize our Python implementation of the sparse Gibbs sampler (with Polya-Gamma Augmentation) for interested researchers.

In addition, we include relevant codes that can replicate every experiment described in our paper.

# Why Our PXL-EM Algorithm?

Here are the main reasons why our proposed Bayesian EM approach can be particularly advantageous for estimating sparse factor loadings from Multidimensional Item Response Theory (MIRT) models:

* **Computational Efficiency**: First, our proposed approach scales much better than MCMC sampling strategies. Unlike the E-step in MCEM or MHRM algorithms which resort to MCMC sampling and cannot be parallelized, we derive the posterior distribution of latent factors explicitly as a type of unified-skew normal distribution inspired by recent result (Durante, 2019). The algorithm is able to sample efficiently from the posterior using minimax tilting methods (Botev,2016). The sampling procedures in the E-step can be parallelized for each individual, and the M-step can also be decomposed into solving $J$ independent penalized probit regressions which can also be estimated efficiently with a standard Glmnet estimation algorithm. Our empirical results show that the algorithm scales well with the numbers of observations, items, and latent dimensions.
* **Latent Dimensionality**: in an unconstrained exploratory MIRT setting, estimating the number of specified factors remains challenging due to the existence of multiple equivalent solutions. In consequence, the number of factors usually has to be specified arbitrarily before fitting the model. We propose to specify an Indian Buffet Process (IBP) prior (Ghahramani and Griffiths, 2015) on the loading matrix, which would allow the algorithm to automatically learn the dimensionality of factors from observed data. Researchers simply need to specify an arbitrarily large truncated level $K^*$ with very little computational cost.
* **Identification & Sparsity**: It is well known that factor analysis suffers from identification issues due to rotational indeterminacy of the loading matrix. In consequence, the estimated factor loading matrix of many exploratory MIRT models tend to be non-sparse and are difficult to interpret. The traditional approach to enhance interpretability is to apply a varimax rotation of the estimated solution after fitting the model. Similar to the approach in (Rockova and George, 2016), we apply the spike-and-slab LASSO (SSL) prior (Rockova, 2018) on the loading matrix, and combine it with an auxiliary rotation matrix induced from parameter expansion to encourage a sparse representation of factor loadings. This approach has the potential to alleviate the inherent identification issues of factor models and produce interpretable loading estimation. Furthermore, this approach can produce exact zeros, which does not need any additional step of thresholding required by MCMC sampling.
* **Robust Initialization**: EM-type algorithms tend to be sensitive to poor initialization and are prone to get stuck in local optimum. We add a factor rotation step after the M-step, which would accelerate EM convergence via parameter expansion. The addition of parameter expansion step (Liu, Rubin, and Wu ,1988) allows the algorithm to explore the parameter space of equivalent observed data likelihoods faster and hence is more robust against poor initialization. One can further adjust the SSL prior to calibrate various initialization methods.   
* **Flexible Loading Structure**:  In contrast to the bifactor MIRT model, which enforces stringent orthogonality conditions among secondary factors, our unconstrained Bayesian estimation approach refrains from imposing such rigid structural assumptions on the factor loading matrix. As a result, our approach frequently identifies more intricate sparse solutions, as evidenced by numerous empirical examples outlined in the paper.


# Repo Directions

* **PXL-EM Algorithm**: Our R implementation of the PXL-EM algorithm can be found in ```project/probit_em_mirt.R```. Some potentially useful functions such as making predictions or computing loglikelihood can be found in ```project/probit_em_util.R```, which is not required to run ```probit_em_mirt.R```. To import the PXL-EM R script, in your R-script editor:
    ```r
    setwd("your_work_dirctory") # set working directory to where the probite_em_mirt.R script is
    source("probit_em_mirt.R")
    source("probit_em_util.R")
    ```
  **Note**: you may find the Python version of the EM algorithm in ```src/bayesian_mirt/em/probit_em.py```. However, that's still under development, **DO NOT USE IT YET**. 
* **Sparse Gibbs Sampler**: Our Python implementation of the sparse Gibbs sampler can be found in ```src/bayesian_mirt/mcmc_sampler/mcmc_mirt.py```. You are welcome to use this script directly to fit our Gibbs Sampler for your applications, or you can also try to download this as a python package in your virtual environment: 
    ```bash
    python setup.py install
    pip install -e .
    ```
    Then in your python script editor:
    ```python
    import bayesian_mirt as bmirt
    ```
* **Experiment Replications**: If you are interested in replicating the experiment results of our PXL-EM algorithm in our paper:
  * **IBP Loading Experiment**: check script ```project/synthetic_data_experiments/s02_probit_em_ibp.R```. The corresponding folder contains scripts to fit the Gibbs sampler, the mhrm algorithm, and the bifactor models.
  * **DESE Experiment**: check script ```project/dese_experiment/s02_fit_dese_data.R```. The corresponding folder contains scripts to fit the Gibbs sampler, the mhrm algorithm, and the bifactor models. Note you may need to run ```s01_process_dese_data.R``` before fitting the PXL-EM algorithm, which preprocesses and samples raw DESE dataset.
  * **QOL Experiment**: check script ```project/qol_experiment/s01_em_qol.R```. The corresponding folder contains scripts to fit the Gibbs sampler, the mhrm algorithm, and the bifactor models.
  
# PXL-EM Algorithm Usage Examples & Advice

We provide a detailed guidance on how to apply our sparse Bayesian MIRT approach, the PXL-EM algorithm, to your particular binary factor analysis applications. To achieve this, we will provide step-by-step instructions in order to replicate the Quality-of-life measurement experiment in our paper. 

This QOL dataset comprises responses from $586$ individuals with chronic mental illnesses, who answered $35$ items related to their life satisfaction across seven distinct subdomains (Family, Finance, Health, Leisure, Living, Safety, and Social). Given the presence of a global life satisfaction item that can serve as the primary factor, Gibbons et al. (Gibbons et al., 2007) showcased the advantages of the bifactor model in a confirmatory analysis. In their study, items associated with each subdomain loaded on separate secondary latent factors as well as the primary quality-of-life factor. 

We will assess whether our unconstrained PXL-EM algorithm can unveil a similar latent structure to that of the bifactor model, even without any external guidance. This particular instance presents a significant challenge, as the algorithm must learn an $8$-dimensional latent space (inclusive of the primary component) from a mere $35$ available items. In the confirmatory analysis of Gibbons, factor loading structures are predefined by experts, here our exploratory analysis received no such information but the item response data.

To better understand these instructions, we invite readers to read the modeling section of our paper, as well as the idea of "Dynamic Posterior Exploration", in which we keep the spike-and-slab LASSO penalty parameters $\lambda_1$ at a reasonable value, and gradually increase the values of $\lambda_0$ for better calibration.

##  Code Example: Replicating QOL Experiment in Our Paper

First set your working directory to where the file ```probit_em_mirt.R``` is. Then load necessary packages. Modify ```data_input_dir```  to the directories where you stored the ```QOL.txt``` data.  Modify ```model_output_dir``` to the directories where you intend to save the PXL-EM model. You can download the QOL dataset from the directory ```data/qol/QOL.txt```

It's worth noting that while the original QOL dataset utilized a seven-point ordinal scale, we transformed it into binary responses by dichotomizing between $4$ and $5$, aligning with our proposed approach's focus on binary response data. Here ```dicho_threshold``` represents the dichotomization threshold to transform ordinal data to the binary data. 
```r
setwd("your workding directory")
source("probit_em_mirt.R")
source("probit_em_util.R")
library(pacman)
p_load("tidyverse", "argparser", "dplyr", "magrittr", "here", "hash")

# parse arguments
arguments <- arg_parser("FIT QOL data") %>% 
  add_argument(
    "--data_input_dir", 
    help= "data_input_dir", 
    default= here("data", "qol")
  ) %>%
  add_argument(
    "--model_output_dir", 
    help= "model_output_dir", 
    default= here("models", "qol")
  ) %>%
  add_argument(
    "--dicho_threshold", 
    help= "how to dichotomize data", 
    default= 4
  ) %>%
  parse_args()

# read data and convert to binary response
lines <- readLines(file.path(arguments[["data_input_dir"]],"QOL.txt"), warn=FALSE)
cleaned_lines <- trimws(lines, "both")
data <- data.frame(cleaned_lines, stringsAsFactors = FALSE)
data <- data %>%
  tidyr::separate(col = cleaned_lines, into = paste0("col", 1:35), sep = 1:34, convert = TRUE)
binary_response <- +(data > arguments[["dicho_threshold"]]) %>% as.matrix() %>% unname()
```
Once we transform the raw QOL data to the binary response data, we will be ready to fit our PXL-EM algorithm. Here:
* ```lambda1```: $\lambda_1$ parameter for the SSL prior. It's usually easier for EM to convergence by specifying a smaller value.
* ```Lambda0_path```: is the $\lambda_0$-path for dynamic posterior exploration. The key is to start with the same value as $\lambda_1$ and gradually increase the $\lambda_0$ values until convergence.
* ```large_k```: serves as an upper bound of our guess on the latent dimensionality (true dimensionality is 8, here we are guessing it's not more than 10, which is reasonable guess for mere 35 items).  
* ```loading_start_large```: random staring positions of the model parameters. ```alphas``` represents the loading matrix, for which it's recommended to start with small positive values, as MIRT models generally have positive loadings. ```intercept``` represents the intercept parameters for the 2PL-MIRT model. ```c_params``` represents the ordered inclusion probabilities for each column of the loading matrix, starting with probabilities a half generally works fine.
* ```dynamic_posterior_exploration``` is the main function fitting our PXL-EM algorithm. Internally, it calls the function ```probit_em_mirt``` function for each value of $\lambda_0$ along the ```lambda0_path``` until convergence.

It took me about 411 seconds to fit our PXL-EM algorithm on my MacBook Pro. 
```r
# fit px-em
lambda1 <- 0.1
lambda0_path <- c(0.1, 0.5, 1, 5, 10, 20, 30, 40)
large_k <- 10
nitems <- 35
set.seed(1)
loading_starts_large <- hash("alphas"= matrix(runif(nitems*large_k, 0, 0.1), nitems, large_k),
                             "intercepts"= runif(nitems , -0.1,0.1), "c_params" = rep(0.5, large_k))
px_em <- dynamic_posterior_exploration(data =binary_response, k = large_k, ibp_alpha = 2, mc_samples =50,
                                       ssl_lambda0_path = lambda0_path, ssl_lambda1 = lambda1, pos_init =TRUE,
                                       max_iterations= 100, epsilon = 0.04, PX = TRUE, varimax = FALSE,
                                       loading_constraints= NULL, start = loading_starts_large,
                                       plot=FALSE, stop_rotation=100, random_state = 1, cores=8)

```
If you want to visualize the estimated factor loading matrix, you can run the function below to plot the loading matrix:
```r
plot(px_em$lambda0_40$alphas, main= "PXL-EM", key=NULL,  axis.row=NULL, axis.column = NULL,  xlab='', ylab='', border=NA) 
```

![](qol_example_plot.png)

If you reorder the column, you will get the same plot in our paper. If you fit a sparse gibbs sampler (see script ```project/qol_experiment/s02_gibbs_qol.py```) or the MHRM algorithm of Cai (see script ```project/qol_experiment/s03_mhrm_qol.R```), you will get loading matrices visualized as below (copied from our paper). The plot highlights the advantage of our PXL-EM algorithm over the other existing exploratory factor analysis strategies, which can produce sparse and interpretable loading estimations that can mirror the oracle bifactor model. Unlike the bifactor model which requires predefined factor loading matrix, we emphasize again our PXL-EM algorithm does not require any external guidance, and can learn the loading structure from the binary item response data. We invite interested readers to check out our paper to learn more about the discussions on this QOL experiment.

![](paper_qol_plot.png)

## Function Documentation

We explain the arguments of ```dynamic_posterior_exploration``` function contained in ```project/probit_em_mirt```.

```r
dynamic_posterior_exploration <- function(data, 
                                          k, 
                                          ibp_alpha, 
                                          mc_samples, 
                                          ssl_lambda0_path, 
                                          ssl_lambda1, 
                                          pos_init =TRUE, 
                                          max_iterations, 
                                          epsilon, 
                                          PX = FALSE, 
                                          varimax = FALSE, 
                                          loading_constraints= NULL, 
                                          start = NULL, 
                                          plot=FALSE, 
                                          stop_rotation=100, 
                                          random_state = 1,
                                          cores=8
                                          )
```

*   ```data```: n by m item response, must be a matrix of zeros and ones.
*   ```k```: trunctaed dimension for the IBP prior. Input your guess of the upper bound of the latent dimension.
*   ```ibp_alpha```: IBP intensity parameter
*   ```mc_samples```: number of monte-carlo sampls to approximate E step. We found setting ```mc_samples``` as $50$ works well for most applications. 
*   ```ssl_lambda0_path```: A vector of $\lambda_0$ parameters of the SSL prior. We recommend starting with $\lambda_0=\lambda_1$, then gradually increase $\lambda_0$. 
*   ```ssl_lambda1```:  $\lambda_1$ parameter for the SSL prior, usually between 0.1 and 1.
*   ```pos_init```: whether to set negative values as zero along the $lambda_0$ path. We recommend setting it to ```TRUE``` as MIRT models typically assume positive loadings.
*   ```max_iterations```: maximum iterations allowed for EM algorithm. The algorithm terminates after ```max_iterations```. 
*   ```epsilon```: stopping threshold. For each value of $\lambda_0$ in ```ssl_lambda0_path```, the algorithm terminates when the infinite norm of the difference of the two consecutive loading matrices is smaller than epsilon. A value around $0.05$ generally works fine. 
*   ```PX```: Boolean parameter, whether to perform parameter Expansion. We recommend always setting it to ```TRUE```.
*   ```varimax```: whether performing varimax rotation. Empirically, we haven't found varimax works well with PXL-EM algorithm so we recommend setting it to ```FALSE```.
*   ```Loading_constraints```: 3-column matrix, first two columns are positions representing row and column indices, last column is fixed value. This argument may be helpful when you feel certain that certain position of the loading matrix should remain a specific value. Our algorithm will make sure not to optimize these entries. For instance, if you want to fix the first component of the first item as one, and the second component of the second item as 0, you can input ```matrix(c(1,1,0,2,2,0), 2, 3, byrow = TRUE)``` 
*   ```start```: parameter initialization, use hash object as illustrated in the example above. We recommend starting with small positive values for the loading matrix (```alphas```), reasonable values for the intercept (```intercepts```) ,and $0.5$ for the ordered inclusion parameter (```c_params```). Example: ```hash("alphas"= matrix(runif(nitems*k, 0, 0.1), nitems, k), "intercepts"= runif(nitems , -0.1,0.1), "c_params" = rep(0.5, k))```
*   ```plot```: Boolean Parameter, whether to plot loading matrix after each iteration, useful for experimentation
*   ```stop_rotation```: switch from PX-EM to EM after "stop_rotation" amount of iterations. You have the option to turn off parameter expansion ```PX``` after certain iterations.
*   ```random_state```: set seed
*   ```cores```: number of cores for parallelization. I usually set it to 8 for my macbook pro. As both the E-step and M-step can be parallelized, you can certainly set it higher for faster iterations when you have enough cores.

Note the ```dynamic_posterior_exploration``` function above involves calibrating the SSL prior along the $\lambda_0$-path. Calibration typically would produce better estimation results, but if you would like to run our PXL-EM algorithm without calibration, you can simply run the function ```probit_em_mirt``` below by applying the PXL-EM algorithm for only one time. The meaning of the required parameters remains the same as above, except for the ```probit_em_mirt``` function, you only need to provide a single value of ```ssl_lambda0```, as opposed of a vector of $\lambda_0$ values.

```r
probit_em_mirt <- function(data, 
                           k, 
                           ibp_alpha, 
                           mc_samples, 
                           ssl_lambda0, 
                           ssl_lambda1, 
                           max_iterations, 
                           epsilon, 
                           PX = TRUE, 
                           varimax = FALSE, 
                           loading_constraints= NULL, 
                           start = NULL, 
                           plot=FALSE, 
                           stop_rotation=100,
                           random_state = 1,
                           cores=8
                           )

```



# Appendix: Notes on Sparse Gibbs Sampler 

While we advocate for our PXL-EM algorithm for estimating sparse MIRT models, if you are interested in our implementation of the Sparse Gibbs Sampler, you can find our it in ```project/src/mcmc_sampler/mcmc_mirt.py```, which is well documented. This efficient gibbs sampler utilizes Polya-Gamma augmentation (Polson, Scott, and Windle, 2013) and the adaptive spike-and-slab prior of Ishwaran and Rao, 2015. It was used as a strong baseline model in our paper and you can find more details and its derivations in the Appendix of our paper. 

Here is the main function to create a sparse Gibbs sampler object:
```python
class MCMC_MIRT:
    def __init__(
            self,
            data: typing.Union[np.ndarray, pd.DataFrame],
            dim: int,
            num_samples: int,
            alpha_prior: str = "normal",
            loading_constraints: typing.Dict[typing.Tuple[int, int], float] = {},
            start: typing.Dict[str, np.ndarray] = {},
            random_state: int = 42,
            auto_zero_out: bool = False
    ):

```
- ```data```: can be a pandas dataframe or a numpy array. Must be binary data.
- ```dim:``` latent dimensionality.
- ```num_samples:``` number of MCMC draws. 
- ```alpha_prior```: the prior for the loading matrix alpha. Default is "normal" without sparsity assumption.  The other two options are "ss" (spike-and-slab prior of George and McCulloch, 1993), and "adaptive-ss" (adaptive spike-and-slab prior of Ishwaran and Rao, 2005). We used "adative-ss" for our paper as the benchmark.
- ```Loading_constraints```: see type hints for required data structure. This is helpful if you want to keep certain elements of the loading matrix fixed: for instance, if you want the element in the second row and the third column fixed as 1, you can specify ```loading_constraints= {"(2, 3)" = 1}```.
- ```start```: You have the option to specify the starting points for ```thetas```(latent factors), ```alphas```(loading matrix), and ```intercepts```. See type hints for required data structure. 
- ```random_state```
- ```auto_zero_out```: we had this idea to zero out very small values of the loading matrix along the sampling process. No good theory or empirical success found so far. So we recommend to keep it as ```False```.

For more examples to fit a gibbs sampler, you may checkout the script ```project/qol_experiment/s02_gibbs_qol.py```, in which we fit a sparse "adaptive-ss" Gibbs sampler using the QOL data.


# Reference

* Chuanhai Liu, Donald B. Rubin, and Ying Nian Wu. Parameter expansion to accelerate
em: The px-em algorithm. Biometrika, 85(4):755–770, 1998.
* Daniele Durante. Conjugate bayes for probit regression via unified skew-normal distribu-
tions. Biometrika, 106(4):765–779, aug 2019.
* Edward I. George and Robert E. McCulloch. Variable selection via gibbs sampling. Journal
of the American Statistical Association, 88(423):881–889, 1993.
* Hemant Ishwaran and J. Sunil Rao. Spike and slab variable selection: Frequentist and
Bayesian strategies. The Annals of Statistics, 33(2):730 – 773, 2005.
* Li Cai. High-dimensional Exploratory Item Factor Analysis by A Metropolis–Hastings
Robbins–Monro Algorithm. Psychometrika, 75(1):33–57, March 2010.
* Nicholas G. Polson, James G. Scott, and Jesse Windle. Bayesian inference for logistic mod-
els using p ́olya–gamma latent variables. Journal of the American Statistical Association,
108(504):1339–1349, 2013.
* Robert D. Gibbons, R. Darrell Bock, Donald Hedeker, David J. Weiss, Eisuke Segawa,
Dulal K. Bhaumik, David J. Kupfer, Ellen Frank, Victoria J. Grochocinski, and Angela
Stover. Full-information item bifactor analysis of graded response data. Applied Psycho-
logical Measurement, 31(1):4–19, 2007.
* Veronika Rockova a and Edward I. George. Fast bayesian factor analysis via automatic
rotations to sparsity. Journal of the American Statistical Association, 111(516):1608–1622,
2016
* Veronika Rockova. Bayesian estimation of sparse signals with a continuous spike-and-slab
prior. The Annals of Statistics, 46(1):401 – 437, 2018.
* Zhehan Jiang and Jonathan L Templin. Gibbs samplers for logistic item response models
via the p ́olya–gamma distribution: A computationally efficient data-augmentation strat-
egy. Psychometrika, 84:358–374, 2018.
* Z. I. Botev. The normal law under linear restrictions: Simulation and estimation via
minimax tilting. Journal of the Royal Statistical Society Series B: Statistical Methodology,
79(1):125–148, feb 2016.
* Zoubin Ghahramani and Thomas Griffiths. Infinite latent feature models and the indian
buffet process. In Y. Weiss, B. Sch ̈olkopf, and J. Platt, editors, Advances in Neural
Information Processing Systems, volume 18. MIT Press, 2005.

