#' AFSSEN_CV
#'
#' For a simulation data, with computing the optimum penalty parameters, It computes important variables and produce smooth estimates of their parameters in a function-on-scalar linear model
#' with sub-Gaussian errors and high-dimensional predictors.
#'
#' @param type_kernel string. three possible choices are implemented.  \code{gaussian},
#' Defualt is \code{"exponential"}.
#' \code{exponential}, \code{sobolev}.
#' @param param_kernel scalar. Value of the characteristic smoothing parameter of the kernel.
#' It is the \eqn{\sigma} parameter of the Gaussian and the Exponential kernel, as introduced
#'  in \link[kernlab]{rbfdot} and \link[kernlab]{laplacedot} functions;
#'  the \eqn{\sigma} parameter of  the Sobolev.
#' Defualt is \code{8}.
#' @param thres scalar. Stopping criteria: beta increment threshold
#' \deqn{
#' || \beta^{T} - \beta^(T-1) ||_{H} < thres
#' }
#' Defualt is \code{0.02}.
#' @param number_non_zeros scalar. Stopping Criteria: Kill switch; number of nonzero predictors
#' Defualt is \code{20}.
#' @param ratio_lambda_H scalar. \eqn{\lambda_{Hmax}/\lambda_{Hmin}}
#' Defualt is \code{0.01}.
#' @param number_lambda_H scalar. Generate the number of log-equally spaced \eqn{\lambda_{H}} in \eqn{[\lambda_{Hmin},\lambda_{Hmax}]}.
#' Defualt is \code{100}.
#' @param num_lambda_H_NONad scalar. Number of \eqn{\lambda_H} in non-adaptive step
#' Defualt is \code{50}.
#' @param lambda_H vector. You have option to insert directly a vector of \eqn{\lambda_H}.
#' Defualt is \code{numeric()}.
#' @param lambda_K vector. Vector of \eqn{\lambda_{K}}.
#'
#' @param early_CV binary. 0 or 1 : applying the \code{early_CV_thres} stopping criteria or not.
#' Defualt is \code{0}.
#' @param early_CV_thres scalar. Stopping Criteria: Breaking point in CV plot.
#' \deqn{
#' |CV(h-1,k)-CV(h,k)| / |CV(h-1,k)| < early_CV_thres
#' }
#' Defualt is \code{0.001}.
#' @param max_ite_nadp scalar. Stopping Criteria: Maximum iteration of coordinate descent algorithm in non-adaptive step
#' Defualt is \code{10}.
#' @param max_ite_adp scalar. Stopping Criteria: Maximum iteration of coordinate descent algorithm in adaptive step
#' Defualt is \code{30}.
#' @param max_ite_final scalar. Stopping Criteria: Maximum iteration of coordinate descent algorithm in final step
#' Defualt is \code{50}.
#' @param target_inc binary. Stopping Criteria: 0 or 1 : if target function is increased, stop
#' Defualt is \code{1}.
#' @param proportion_training_set scalar. value in (0,1), the
#' proportion for the training set for the Cross Validation in non-adaptive step
#' Defualt is \code{0.75}.
#' @param fold_ad scalar. Number of fold for using CV in adaptive steps to find optimum \eqn{\lambda_{H}} and \eqn{\lambda_{K}} and then the coefficients estimation.
#' Defualt is \code{10}.
#'
#' @return  list containing:
#'  \itemize{
#'  \item \code{beta : } matrix. final estimation of coefficients.
#'  \item \code{beta_no_adaptive : }  matrix. estimation of coefficients after non-adaptive step.
#'  \item \code{predictors : } vector. final significant predictors.
#'  \item \code{predictors_no_adaptive : } vector. significant predictors after non-adaptive step.
#'  \item \code{lambda_H_opt : } scalar. optimum \eqn{\lambda_{H}}
#'  \item \code{lambda_K_opt : } scalar. optimum \eqn{\lambda_{K}}
#'  \item \code{pred_error : } scalar. prediction error
#'  \item \code{pred_error_der : } scalar. prediction error derivative.
#'  }
#'
#' @export
#' @examples
#' \dontrun{
#' data(simulation)
#' data(SobolevKernel)
#' time <- proc.time()
#' FLAME_estimation <- FLAME()
#' duration <- proc.time()-time
#' duration
#' names(FLAME_estimation)
#' }
#'

AFSSEN_CV <- function(type_kernel="exponential",param_kernel=8,
                     thres, number_non_zeros,
                     ratio_lambda_H, number_lambda_H, num_lambda_H_NONad,
                     lambda_H, lambda_K,
                     early_CV, early_CV_thres, max_ite_nadp, max_ite_adp, max_ite_final, target_inc,
                     proportion_training_set=0.75, verbose=FALSE, fold_ad=10,
                     N=500, I=1000, I0=10,
                     nu_beta=2.5, range=1/4, variance=1, # for the rough beta in Matern
                     nu_eps = 1.5, range_eps=1/4 # for epsilon in Matern
                     )
{

############################################################################
#
# Part 0: Defining the function for numerically estimate norm beta in H
#
############################################################################

    estimation_norm_R <- function(lambda_H, lambda_K, omega, B, tau)
    {

        ## @A: the function we want to optimize in AFSSEN is:
        ## 1=\sum_{j=1}^{\infty} \dfrac{tau_j < \beta_j , v_j >^2 }{ tau_j \norm{\beta_j}_K + \lambda omega_j}^2


        ## @A: nloptr is an R package for nonlinear optimization.

        ## ?@A: why do you optimiza |1-1/sum| instead of |1-sum|?

        optimization_function <- function(x, lambda_H, lambda_K, omega, B, tau)
        {
            numerat <- (B^2) * (tau^2)
            denom <- ((tau+lambda_K)*x + lambda_H*omega*tau)^2
       #    tot <- 1 - 1 / (sum(numerat/denom))
            tot <- 1 - sum(numerat/denom)
            return(abs(tot))
        }

        ## ?@A: Are these initial values are fixed or we need to change them in different situation?
        opts= list("algorithm"= "NLOPT_LN_COBYLA", "xtol_rel"=1.0e-16)


        ## ?@A: why it ended up with semicolon (;)?
        ott_model <- nloptr(0, optimization_function, opts=opts, lb=0,
                            lambda_H=lambda_H,lambda_K=lambda_K, omega=omega, B=B, tau=tau);

        return(ott_model$solution)

    }


    Proc_Time <- function()
    {
        r=as.numeric(proc.time()[3])
        return(r)
    }

    function_computation <- estimation_norm_R
    function_time <- Proc_Time



############################################################################
#
# Part 1: Data generation
#
############################################################################



    m <- 50 # total number of points
    T_domain <- seq(0, 1, length = m) # time points, length = m
    M_integ <- length(T_domain)/diff(range(T_domain))

    # # # # # #######
    # # #
    # # defintion of the design matrix X, in this specific case the
    # # covariance matrix C is the identity matrix
    mu_x <- rep(0, I)
    C <- diag(I)
    X <- mvrnorm(n=N, mu=mu_x, Sigma=C) # X is a N*I matrix
    X <- scale(X) # normalization

    # # # # # #######
    # # # #
    # # defintion of the coefficients

    hyp <- c(log(range), log(variance)/2) # set of parameters for the # Matern Covariance operator of beta
    mu_beta <- rep(0,m) # mean of the beta
    Sig_beta <- covMaterniso(nu_beta, rho = range, sigma = sqrt(variance) , T_domain)
    beta <- mvrnorm(mu=mu_beta, Sigma=Sig_beta, n=I0) # beta is a I0*m matrix. Each row is made by a matern process.

    # # # # # # I0 significant coefficients
    # # # # # # defintion of the random errors

    mu_eps <- rep(0, m)
    Sig_eps <- covMaterniso(nu_eps, rho = range_eps, sigma = sqrt(variance), T_domain)
    eps <- mvrnorm(mu=mu_eps, Sigma=Sig_eps, n=N) # eps is a N*m matrix

    # # # # # random errors
    # # #
    I_X <- sort(sample(1:I, I0)) # index of the I0 significant predictors

    Y_true <- X[,I_X]%*%beta
    Y_full <- X[,I_X]%*%beta + eps # Y_n observations (pointwise evaluation)



    ### Generating Kernel

    M <- 50

    kernel_here <- generation_kernel(type = type_kernel,parameter = param_kernel,
                                     domain = T_domain, thres = 0.99, return.derivatives = TRUE)

    eigenval <- kernel_here$eigenval
    eigenvect <- kernel_here$eigenvect
    derivatives <- kernel_here$derivatives


    # # # # preojection on the kernel basis of y and beta
    Y_matrix <- projection_basis(Y_full, eigenvect, M_integ)
    XB <- projection_basis(Y_true, eigenvect, M_integ)


    B_true <- projection_basis(beta, eigenvect, M_integ)
    matrix_beta_true_full <- matrix(0, dim(B_true)[1], I)
    matrix_beta_true_full[,I_X] <- B_true

    # # #derivatives
    B_true_der_time <- t(kernel_here$derivatives %*% B_true)
    Y_true_der_time <- X[,I_X]%*%B_true_der_time
    XBd <- projection_basis(Y_true_der_time, eigenvect[-dim(eigenvect)[1],], M_integ-1)




    Y <- Y_matrix




############################################################################
#
# Part 2: RUN non-adaptive version to estimte weights for adaptive step
#
############################################################################

#######################################################
#
# Part 2.1: estimation_first
#
#######################################################

    weights = rep(1, dim(X)[2])

    # First step: all the weights are set to 1
    # first run to identify the possible values of lambda
    # (from lambda_max which makes all the predictors set to 0, to the the lambda value which
    # number_non_zeros predicotors are different from 0)

    # num_lambda_H_NONad = number_lambda_H  # NONad= Non adaptive

    #num_lambda_NONad <- round(number_lambda/5) # for this first step (the most
    # computationally expensive) we can reduce the number of lambda examined.
    # Then (in the Second: Adaptive Step)
    # we will consider the whole set of lambdas to define the final estimation






#######################
#######################
#######################     Non-adaptive
#######################
#######################







        if (verbose) print('Non-adaptive step')

        subset <- c(rep(2, ceiling(dim(X)[1]*proportion_training_set)),
                    rep(1, ceiling(dim(X)[1]*(1-proportion_training_set))))

        # proportion_training_set is the percentage of the data to use to
        # fit the model (2), the remaining part to compute the prediction
        # error:left_out (1)

        random_groups <- subset[sample(1:dim(X)[1])]



        # training and test set
        i <- 1

        left_out <- which(random_groups==i)

        # fitting of the model with the proportion_training_set% of the data and
        # computation of the CV error

        #sourceCpp("FLAME_functions_cpp.cpp")

        eigenvect_for_der=eigenvect[-dim(eigenvect)[1],]


        estimation_first_CV <- definition_beta_CV(X[-left_out,], Y[, -left_out],
                                 X[left_out,], Y[, left_out],
                                 eigenval, weights,
                                 function_computation, function_time,
                                 thres, number_non_zeros,
                                 lambda_H,lambda_K,
                                 ratio_lambda_H, num_lambda_H_NONad,
                                 early_CV, early_CV_thres, max_ite_nadp, target_inc,
                                 verbose, base::sum)
                                 #projection_basis, M_integ-1,
                                 #derivatives, eigenvect_for_der, XBd[,left_out])




        ## @A: estimation_first_CV$error is the \norm{Y_pred - Y_test}_H for different lambda.
        ## @A: lambda_first_selected is not the optimum lambda. It is the index of the lambda.

        ## @A: estimation_first_CV$error is a matrix (#lambda_H)*(#lambda_K)
        lambda_H_first_selected <- which(estimation_first_CV$error==min(estimation_first_CV$error[estimation_first_CV$error>0]) , arr.ind = TRUE)[1] # optimum lambda_H
        lambda_K_first_selected <- which(estimation_first_CV$error==min(estimation_first_CV$error[estimation_first_CV$error>0]) , arr.ind = TRUE)[2] # optimum lambda_K


#######################################################
#
# Part 2.3: estimation_first_definite, weights_new and
#
#######################################################


        if (length(estimation_first_CV$Pred_opt)==0)
        {
            print('No significant predictors indentified.')
            result <- list(beta=NULL,
                           beta_no_adaptive=NULL,
                           predictors=NULL,
                           predictors_no_adaptive=NULL)
            return(result)
        }else
            {
                if(length(estimation_first_CV$Pred_opt)==1)
                {

                    ## @A: In the following matrices, we have just one column because the length(estimation_first_definite$Pred)==1
                    beta_selected <- matrix(estimation_first_CV$Beta_opt[,estimation_first_CV$Pred_opt],
                                            length(estimation_first_CV$Beta_opt[,estimation_first_CV$Pred_opt]), 1)
                    X2_data <- matrix(X[,estimation_first_CV$Pred_opt], length(X[,estimation_first_CV$Pred_opt]),1)

                }

            else{
                beta_selected <- estimation_first_CV$Beta_opt[,estimation_first_CV$Pred_opt]
                X2_data <- X[,estimation_first_CV$Pred_opt]
            }
            }



            # defintion of the weights
            weights_new <- 1/norm_matrix_H(beta_selected)

            ## @A: For adaptive step, we just consider the non zero predictors and delete the zero ones.
            ## you can see the length of weights_new can be less than the initial weights.

############################################################################
#
# Part 3: RUN adaptive step
#
############################################################################


#######################################################
#
# Part 3.1: estimation_second
#
#######################################################


            if (verbose) print('Adaptive step')


            flds <- createFolds(1:dim(X)[1], k = fold_ad, list = TRUE, returnTrain = FALSE)

            SUM_MAT_CV = matrix(0,number_lambda_H,length(lambda_K))
            SUM_MAT_der = matrix(0,number_lambda_H,length(lambda_K))
            #error_der = rep(0,fold_ad)

      for( f in 1:fold_ad){


      left_out <- flds[[f]]



      estimation_second_CV <- definition_beta_CV(as.matrix(X2_data[-left_out,]), Y[, -left_out],
                                        as.matrix(X2_data[left_out,]), Y[, left_out],
                                        eigenval, weights_new,
                                        function_computation, function_time,
                                        thres, number_non_zeros,
                                        lambda_H, lambda_K,
                                        ratio_lambda_H, number_lambda_H,
                                        early_CV, early_CV_thres, max_ite_adp,target_inc,
                                        verbose, base::sum)
                                        #projection_basis, M_integ-1,
                                        #derivatives, eigenvect_for_der, XBd[,left_out])


      A <- estimation_second_CV$error

      A [A ==0] <- Inf

      SUM_MAT_CV = SUM_MAT_CV + A

      # B <- estimation_second_CV$error_der
      #
      # B [B ==0] <- Inf
      #
      # SUM_MAT_der = SUM_MAT_der + B


      # # # #derivatives
      # B_pred_der_time <- t(derivatives %*%estimation_second_CV$Beta_opt[,estimation_second_CV$Pred_opt])
      # Y_true_der_time <- X[,estimation_second_CV$Pred_opt]%*%B_pred_der_time
      # Y_pred_der <- projection_basis(Y_true_der_time, eigenvect[-dim(eigenvect)[1],], M_integ)
      #
      # diff_der= Y_true_der - Y_pred_der
      #
      # error_der[f]=sqrt(sum(apply(diff_der,2,function(x) sum(x^2))))

      print(paste("fold",f,"   ,"))
            }


            Mean_Error_CV=as.matrix(SUM_MAT_CV/fold_ad)
            # Mean_Error_der=as.matrix(SUM_MAT_der/fold_ad)

            ind_min_cv = as.vector(which(Mean_Error_CV==min(Mean_Error_CV[which(Mean_Error_CV>0)]) , arr.ind = TRUE)[1,])


            min_pred_CV=min(Mean_Error_CV)

            # if(length(lambda_K)==1)
            # corres_pred_der=Mean_Error_der[ind_min_cv[1]]
            # else
            # corres_pred_der=Mean_Error_der[ind_min_cv]
            #
            # min_pred_der=min(Mean_Error_der)

            lambda_K_opt=lambda_K[ind_min_cv[2]]
            lambda_H_opt=estimation_second_CV$lambda_H[ind_min_cv[1]]



###################################
#
# Run for the whole data with best best lambda_H and lambda_K
#
###################################


            if (verbose) print('Adaptive step for optimum parameters')

            estimation_second <- definition_beta(as.matrix(X2_data), Y,
                                        eigenval, weights_new,
                                        function_computation, function_time,
                                        thres, number_non_zeros,
                                        lambda_H_opt, lambda_K_opt,
                                        ratio_lambda_H, number_lambda_H,
                                        verbose,base::sum)





            # if (verbose) print('Adaptive step')
            #
            # subset <- c(rep(2, ceiling(dim(X)[1]*proportion_training_set)),
            #             rep(1, ceiling(dim(X)[1]*(1-proportion_training_set))))
            #
            # # proportion_training_set is the percentage of the data to use to
            # # fit the model (2), the remaining part to compute the prediction
            # # error:left_out (1)
            #
            # random_groups <- subset[sample(1:dim(X)[1])]
            #
            #
            #
            # # training and test set
            # i <- 1
            #
            # left_out <- which(random_groups==i)
            #
            #
            # estimation_second_CV <- definition_beta_CV(as.matrix(X2_data[-left_out,]), Y[, -left_out],
            #                             as.matrix(X2_data[left_out,]), Y[, left_out],
            #                             eigenval, weights_new,
            #                             function_computation, function_time,
            #                             thres, number_non_zeros,
            #                             lambda_H, lambda_K,
            #                             ratio_lambda_H, number_lambda_H,
            #                             early_CV, early_CV_thres, max_ite_adp,target_inc,
            #                             verbose, base::sum)



#######################################################
#
# Part 3.3:  estimation_second_definite and finding lamda which implies minimum CV
#
#######################################################


#######################################################
#
# Part 3.3:  The output values
#
#######################################################

            predictors_2=estimation_second$Pred # final set of
            # predictors different from 0 (among the ones isolated in the
            # First: Non-Adaptive step)
            beta_2=estimation_second$Beta[,predictors_2] # estimated betas,
            # it is the matrix of the coefficients of the basis expansion
            # of the betas (with respect to the basis defined by the eigenvectors
            # of the kernel)

            ## @A: Pay attention: predictors_2 are the labels of nonzro predictors among the
            ## ones isolated in the First: Non-Adaptive step. So for finding the actual labels among the
            ## whole predictors we need to define "predictor_def" as follows

            predictor_def <- estimation_first_CV$Pred_opt[predictors_2] # index of
            # the non-zero predictors computed in the estimation

            beta_def <- matrix(0, dim(estimation_first_CV$Beta_opt)[1],
                               dim(estimation_first_CV$Beta_opt)[2])

            ## @A: "beta_def" is a J*I matrix which the nonzero ones estimated by estimation_second and
            ## rest of the columns are zero.

            lambda_H=estimation_second_CV$lambda_H
            lambda_K=estimation_second_CV$lambda_K






            beta_def[, predictor_def] <- beta_2 # final matrix of the estimated betas



            XBhat=(beta_2)%*%t(X[,predictor_def])


            # # #derivatives
            Bhat_true_der_time <- t(kernel_here$derivatives %*% beta_2)
            Yhat_true_der_time <- X[,predictor_def]%*%Bhat_true_der_time
            XBhatd <- projection_basis(Yhat_true_der_time, eigenvect[-dim(eigenvect)[1],], M_integ-1)


            pred_error=sum(norm_matrix_H(XB-XBhat))
            pred_error_der=sum(norm_matrix_H(XBd-XBhatd))

            # (still as coefficients of the kernel basis)
            result <- list(beta=beta_def,
                           beta_no_adaptive=estimation_first_CV$Beta_opt,
                           predictors=predictor_def,
                           predictors_no_adaptive=estimation_first_CV$Pred_opt,
                           true_pred=I_X,
                           Mean_Error_CV=Mean_Error_CV,
                           min_pred_CV=min_pred_CV,
                           index_lambdas=ind_min_cv,
                           lambda_H=lambda_H,
                           lambda_K=lambda_K,
                           lambda_K_opt=estimation_second_CV$lambda_k_opt,
                           lambda_H_opt=estimation_second_CV$lambda_H_opt,
                           L_matrix_nonad=estimation_first_CV$L_matrix,
                           L_matrix_ad=estimation_second_CV$L_matrix,
                           Mat_target_inc=estimation_second_CV$Mat_target_inc,
                           vec_number_lambda_H_computed_adp=estimation_second_CV$vec_number_lambda_H_computed,
                           vec_break_in_number_pred_adp=estimation_second_CV$vec_break_in_number_pred,
                           vec_break_in_early_cv_adp=estimation_second_CV$vec_break_in_early_cv,
                           vec_norm_H_diff_beta_inc_nadp=estimation_first_CV$vec_norm_H_diff_beta_inc,
                           vec_norm_H_diff_beta_inc_adp=estimation_second_CV$vec_norm_H_diff_beta_inc,
                           X=X,
                           Y=Y,
                           B_true=B_true,
                           beta_time=beta,
                           XB=XB,
                           XBd=XBd,
                           XBhat=XBhat,
                           XBhatd=XBhatd,
                           pred_error=pred_error,
                           pred_error_der=pred_error_der,
                           kernel_vec_der=derivatives,
                           kernel_vec=eigenvect,
                           kernel_val=eigenval,
                           T_domain=T_domain)

    return(result)

      }






