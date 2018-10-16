#' Finding the neighborhood of best lambda_H and lambda_K
#'
#' This function computes the neighborhood of the best lambda_K and lambda_H.
#'
#' @param AFFSEN_CV. output of the AFSSEN_CV function
#' @param diff_between_K. the threshold for finding the optimum K
#' @param diff_within_H. the threshold for finding the optimum H
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{lambda_K_opt}. the optimum lambda_K
#'  \item \code{lambda_H_opt}. the optimum lambda_H
#'
#' @export
#'
#' @useDynLib flm, .registration = TRUE
#'

opt_lambdas<- function(AFSSEN_CV, diff_between_K , diff_within_H)
{


    dim_H=dim(AFSSEN_CV$CV_error)[1]
    dim_K=dim(AFSSEN_CV$CV_error)[2]

##########
##########

    opt=which(AFSSEN_CV$CV_error==min(AFSSEN_CV$CV_error),arr.ind = TRUE)

if((opt[1]!=dim_H)||(opt[2]!=dim_K)){
    print("optimum is not in the boundry.")
    result=c(AFSSEN_CV$lambda_H[opt[1]] , AFSSEN_CV$lambda_K[opt[2]])
    names(result)=c("lambda_H_opt","lambda_K_opt")
    return(result)
} else{

##########
##########

    diff_K=rep(NA,dim(AFSSEN_CV$CV_error)[2]-1)

    for(i in 1:(dim(AFSSEN_CV$CV_error)[2]-1)){
        diff_K[i]=mean((AFSSEN_CV$CV_error[,i+1]-AFSSEN_CV$CV_error[,i])^2)
    }

    if(min(diff_K) < diff_between_K){
        ind_K=min(which(diff_K < diff_between_K))
        lambda_K_opt=AFSSEN_CV$lambda_K[ind_K]
    }else{
        print("You should choose more value for lambda_K")
        ind_K=dim(AFSSEN_CV$CV_error)[2]
        lambda_K_opt=AFSSEN_CV$lambda_K[ind_K] # return the last lambda_K
    }

##########
##########

    opt_H=which(AFSSEN_CV$CV_error[,ind_K]==min(AFSSEN_CV$CV_error[,ind_K]),arr.ind = TRUE)

    if(opt_H!=dim_H){
        print("optimum lambda_H is not in the boundry.")
        result=c(AFSSEN_CV$lambda_H[opt_H] , lambda_K_opt)
        names(result)=c("lambda_H_opt","lambda_K_opt")
        return(result)
    }else{

    diff_H=rep(NA,dim(AFSSEN_CV$CV_error)[1]-1)


    for(i in 1:(dim(AFSSEN_CV$CV_error)[1]-1)){
        diff_H[i]=abs((AFSSEN_CV$CV_error[i+1,ind_K]-AFSSEN_CV$CV_error[i,ind_K])/(AFSSEN_CV$CV_error[i+1,ind_K]-AFSSEN_CV$CV_error[dim_H,ind_K]))
    }

    if(min(diff_H[diff_H>0]) < diff_within_H){   # min nonzero diff_H
        ind_H=min(which( diff_H[diff_H>0] < diff_within_H ))
        lambda_H_opt=AFSSEN_CV$lambda_H[ind_H]
    }else{
        print("You should choose more value for lambda_H")

        lambda_H_opt=tail(AFSSEN_CV$lambda_H,1) # return the last lambda_K
    }


##########
##########

    result=c(lambda_K_opt, lambda_H_opt)
    names(result)=c("lambda_K_opt","lambda_H_opt")

    return(result)
    }
    }


    }

