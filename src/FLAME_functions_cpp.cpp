#include <math.h>
#include <cstdio>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <functional>
#include <time.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

using namespace Rcpp;
using namespace arma;
using namespace std;



class linear_model {



  /*###########################################################################################
#
# @A Part0: Defining private values
#
###########################################################################################*/

  // @A: A class (user_defined types) is used to specify the form of an object and it combines
  // data representation and methods for manipulating that data into one neat package.

private:

  // @A: arma (Armadillo) is a high quality linear algebra library (matrix maths) for the C++
  // arma::mat the main matrix object in arma is mat
  // arma:vec  define a vector
  arma::mat Y, X, B, E, B_ls;
  unsigned int num_pred;
  unsigned int num_data;
  unsigned int num_zeros_pred;
  double lambda_H_max;

  /*###########################################################################################
#
# @A Part1: Defining some norm_K and norm_H functions over B and Y
#
###########################################################################################*/

public:

  // c++ function which calls the R function to perform the optimization and
  // compute the norm of beta_j


  // Brief introduction on Rcpp Package:
  // @A: Rcpp package: provides matching C++ classes for a large number of basic R data types.
  // The mapping of data types works in both directions. Transfer from R to C++, and from C++ to R
  // through the templated functions Rcpp::as (for conversion of objects from R to C++) and Rcpp::wrap (for conversion from C++ to R)
  // Hence, a package author can keep his data in normal R data structures without having to worry
  // about translation or transfering to C++. At the same time, the data structures can be accessed as easily at the C++ level, and used in the normal manner.
  // R data types (SEXP : S-expression pointer) are matched to C++ objects in a class hierarchy.
  //  For example, numeric vectors are represented as instances of the Rcpp::NumericVector class,
  // environments are represented as instances of Rcpp::Environment,
  // functions are represented as Rcpp::Function, etc ...
  // Rcpp::wrap  The underlying C++ library also offers the Rcpp::wrap function which is a templated function
  // that transforms an arbitrary object into a SEXP (R type). We should  wrap them when they need to be returned to R.
  // Rcpp::as    The reverse conversion (from R into C++) is performed by the Rcpp::as
  // // conversion from R to C++:
  // template <typename T> T as( SEXP x) ;
  // conversion from C++ to R:
  // template <typename T> SEXP wrap(const T& object) ;


  // @A: omega is the weights in AFSSEN function.
  // lambda_H is the sparsity penalty paramter (\lambda_{2}) in the AFSSEN function
  // lambda_K is the smoothing penalty paramter (\lambda_{1}) in the AFSSEN function
  // computation_norm_cpp is used to find the \norm{\beta_j}_K numerically from the equation.
  // and don't forget the f function is gonna be come from R.



  // @@A: For adding the effects of the L operator, we should add another vector "eta" for its eigenvalues.
  // Here we consider L=identity

  double computation_norm_cpp( double &lambda_H,double &lambda_K , double &omega, arma::vec &B_j,
                               arma::vec &tau,  Rcpp::Function f)
  {
    return(Rcpp::as<double>(f(lambda_H,lambda_K,omega, B_j, tau)));
  }


  double Time(Rcpp::Function T)
  {
    return(Rcpp::as<double>(T()));
  }


  // @A: K is a kernel with k eigenvalues (tau) and k eigenfunc (V). So we have tau=[tau_1,...,tau_k] and V=[v_1,...,v_k]
  // B' is coefficient of the model which are func so we can fix the V as basis and write B'=c_1v_1+...+c_kv_k
  // So B in this code is a k vector including its basis coefficients B=[c_1,...,c_k]
  // and if there are p predictors, there will be p coefficients B'1,...B'p
  // So B_here is a k*p matrix. Each column including their basis coefficients.
  // X_here is a N*p matrix
  // Y_here is a k*N. Each column is the basis coefficients of Y_n. Because when Y_n=X_1B'1+...X_pB'p, \alpha_ni=<Y_n,v_i>=X_1c_1i+...X_pc_pi which is come from B*t(X)
  // p is the I in the orignial model (p=I)

  // norm in K of K(x), given the components in the kernel basis (B) and the
  // eigenvalues of the kernel (tau)


  // @A: norm_K_Kx computes \norm{K(B)}_K ; B is a coefficient in model
  // \norm{k(B)}^2_K = \sum_{i=1}^{k} \dfrac{<K(B),v_i>^2}{\tau_i}=\sum_{i=1}^{k} \dfrac{<B,K(v_i)>^2}{\tau_i}=\sum_{i=1}^{k} tau_i c_i^2
  // when we compute norm_K we should have both eigenval and eigenfunc but in here tau is enough. Becauase the effects of eigenfunc of kernel is in the coefficient in matrices B_here and Y_here.
  //  double norm_K_Kx(arma::vec &B, arma::vec &tau)
  //  {
  //    return(sqrt(sum(B % B % tau)));
  //  }

  // K norm of the matrix B_here (square-root of the sum of the squares of the
  // K norms of the columns of B_here). tau is the vector of the eigenvalues of
  // the kernel

  // @A: norm_K_matrix computes \sqrt{ \norm{B_1}^2_K + ... + \norm{B_p}^2_K } ; B_i is a coefficient in model
  double norm_K_matrix(arma::mat &B_here, arma::vec tau)
  {
    double norm_matrix=0;
    for (unsigned int i=0; i<B_here.n_cols; i++)
    {
      norm_matrix += sum(B_here.col(i)%B_here.col(i)/tau);
    }
    return (sqrt(norm_matrix));
  }

  // square of the H norm of the matrix Y_here (sum of the squares of the
  // H norms of the columns of Y_here)

  // @A: square_norm_H_matrix computes { \norm{Y_1}^2_H + ... + \norm{Y_N}^2_H }
  // \norm{Y_1}^2_H= \alpha_11^2+...+\alpha_1k^2 (elements in first column of Y_here).
  double square_norm_H_matrix(arma::mat &Y_here)
  {
    double norm_matrix=0;
    for (unsigned int i=0; i<Y_here.n_cols; i++)
    {
      norm_matrix += sum(Y_here.col(i)%Y_here.col(i));
    }
    return (norm_matrix);
  }


  double norm_H_matrix(arma::mat &Y_here)
  {
    double norm_matrix=0;
    for (unsigned int i=0; i<Y_here.n_cols; i++)
    {
      norm_matrix += sum(Y_here.col(i)%Y_here.col(i));
    }
    return (sqrt(norm_matrix));
  }



  // @A: norm_H_vec computes \norm{B}_H=\sqrt{b_1^2+...+b_k^2} where b_i is the projection of B on basis.
  double norm_H_vec(arma::vec &B)
  {

    return(sqrt(sum(B%B)));
  }


  // square of the K norm of the matrix Y_here (sum of the squares of the
  // K norms of the columns of Y_here). tau is the vector of the eigenvalues of
  // the kernel

  // @A: square_norm_K_matrix { \norm{Y_1}^2_K + ... + \norm{Y_N}^2_K } ;
  // \norm{Y_1}^2_K = \sum_{i=1}^{k} \dfrac{\alpha_1i^2}{\tau_i}
  double square_norm_K_matrix(arma::mat &Y_here, arma::vec &tau)
  {
    double norm_matrix=0;
    for (unsigned int i=0; i<Y_here.n_cols; i++)
    {
      norm_matrix += sum(Y_here.col(i)%Y_here.col(i)/tau);
    }
    return (norm_matrix);
  }

  double target_function(arma::mat B, arma::vec tau, double lambda_K, double lambda_H, arma::vec omega_vect){
    double out=0;
    arma::vec B_ite;
    arma::mat Y_BX;
    for (unsigned int s=0; s<num_pred; s++){
      B_ite=B.col(s);
      out=out+omega_vect(s)*norm_H_vec(B_ite);
    }
    Y_BX=Y - B*X.t();
    double result=(1/(2*num_data))*square_norm_H_matrix(Y_BX)+(lambda_K/2)*square_norm_K_matrix(B,tau)+lambda_H*out;
    return(result);
  }


  /*###########################################################################################
#
# @A Part2: Defining BIC
#
###########################################################################################*/

  // computation of the BIC given as input
  // - Y_here: kxN matrix; the matrix of the true values of the observations in the kernel basis
  // - X_here: Nxp matrix; predictors' matrix
  // - B_here: kxp matrix; estimated coefficients
  // - p: int; number of non-zero predictors identified
  // - N: int; number of observations
  // - tau: k vector; eigenvalues of the kernel
  double computation_BIC (arma::mat &Y_here, arma::mat &X_here, arma::mat &B_here,
                          unsigned int p, unsigned int N, arma::vec &tau)
  {
    arma::mat difference = Y_here - B_here*X_here.t();
    // we can compute the RSS both with the H and K norm.
    // double error_H = square_norm_H_matrix(difference); //H norm


    // ?@A: shouldn't it be error_H=square_norm_K_matrix(difference)? Because the Y may not exist in K.
    // double error_K = square_norm_K_matrix(difference, tau); //K norm
    double error_H = square_norm_H_matrix(difference); //H norm

    //std::cout<<error_H<<N<<p<<std::endl;
    double BIC_computed =  N*log(error_H/N) + log(N)*p;
    return BIC_computed;

  }

  /*###########################################################################################
#
# @A Part3: Defining B_ls .... Initializing the private value of the linear_model class.
#
###########################################################################################*/

  // constructor of the linear_model element given the matrix of the predictors
  // X_data (Nxp) and the matrix of the observations Y_data (kxN)

  // ?@A: are these two following commands for initialization private values? when do we usually use it?

  // @A: B_ls is actually the B^{v} in FLAME. B_i^{v}=\dfrac{\sum_{n=1}^{N} X_n,i E_n}{N} where E_n=Y_n-\sum_{j != i} X_n,j B_j
  // Since we can write (AB)_ij=\sum_{k} A_{ik} B_{kj} so B_ij^{v}=\dfrac{\sum_{n=1}^{N} X_n,i E_jn}{N}=\dfrac{\sum_{n=1}^{N}E_jn X_n,i}{N}
  // So we can write B_ls= E*X/N

  // .zeros makes all elements zero.

  // @A: We did not clarify the dimentsion of B before. when we say
  // B.zeros(B_ls.n_rows, B_ls.n_cols), it means B is a B_ls.n_rows*B_ls.n_cols zero matrix.


  linear_model(const arma::mat &X_data, const arma::mat &Y_data){
    Y=Y_data;
    X=X_data;
    E=Y; // inital estimation of the error equal to the data, since beta is initalized with 0
    num_data = X.n_rows;
    B_ls = E * X / num_data; // first LS estimation of beta
    num_pred= X.n_cols;
    num_zeros_pred=0;
    B.zeros(B_ls.n_rows, B_ls.n_cols);
  }

  // constructor of the linear_model element given the matrix of the predictors
  // X_data (Nxp), the matrix of the observations Y_data (kxN) and the first
  // estimation of the coefficients Beta (kxp)


  // ?@A: when B=0 we can define E=Y and say B_ls=E*X/N, but when B=Beta we know
  // B_ls=E*X/N where E_n=Y_n-\sum_{j != i} X_n,j B_j. In the following formula, we defined E=Y - B*X.t()
  // Although I can see in future, you will add B*X_j to E_glob and I guess it handle mu issue here.
  // I mean in this part, B_ls is not B^{v} in AFSSEN function!

  // @A: I guess, we do not use this linear model with Beta, in future!
  linear_model(const arma::mat &X_data, const arma::mat &Y_data, const arma::mat &Beta){
    Y=Y_data;
    X=X_data;
    num_data = X.n_rows;
    B = Beta;
    E= Y - B*X.t();
    B_ls = E * X / num_data;
    num_pred= X.n_cols;
    num_zeros_pred=0;
  }

  // @A: "cout" means "standard output stream" and "<<" meanning "write to"
  //        Rcpp:Rcout works in R as cout in C++


  // output function to print the number of the predictors and the predictor matrix
  void predictors()
  {
    Rcpp::Rcout<<"matrix of predictors" <<X <<"number of predictors "<< num_pred;
  }

  //make the betas 0

  // @A: .zers() : set all elements to zero
  void set_B_zeros()
  {
    B.zeros();
  }

  /*###########################################################################################
#
# @A Part4: definition for choosing lambda_H vector
#
###########################################################################################*/

  // @A: definition of the vector of all the possible values of lambda_H,
  // given the number of lambda_Hs (num_lambda_H) and the ratio between the maximum
  // and the minimum (ratio_lambda_H).
  // in particular the vector is started from lambda_H_max
  // (lambda_H_max = the minimum value of lambda_H which all the predictors are guaranteed to be 0)
  // to lambda_H_min = lambda_H_max*raito_lambda_H
  // not equispaced vector

  // @A: By default we defined them as equispaced in log
  // omega is the vector of weights in AFSSEN function
  arma::vec definition_lambda_H(double &ratio_lambda_H, int &num_lambda_H, arma::vec &omega, arma::vec &tau)
  {

    // @A: lambda_H_max_vec is just an auxiliary vector to help us find lambda_H_max
    arma::vec lambda_H_max_vec(num_pred);  // p*1 vector
    arma::vec lambda_H(num_lambda_H);
    for (unsigned int j=0; j<num_pred; j++)
    {
      arma::vec B_temp=B_ls.col(j);
      lambda_H_max_vec(j) = norm_H_vec(B_temp)/omega(j);

      // ?@A: what is the difference between Rcpp::Rcout and std::cout?

      //std::cout<<"j="<<j<<" lambda "<<lambda_max_vec(j)<<std::endl;
    }


    // @A: lambda_max was defined as private values

    // ?@A: Isn't the lambda_max = 2*max(lambda_max_vec) too high? Can't it be just lambda_max =max(lambda_max_vec)?
    // lambda_H_max = 2*max(lambda_H_max_vec); //CAMPPPPP
    lambda_H_max = max(lambda_H_max_vec); //CAMPPPPP

    // std::cout<<"maximum lambda_H "<<lambda_H_max<<std::endl;
    double lambda_H_min= lambda_H_max*ratio_lambda_H;
    // std::cout<<"defintion of a set of lambda_H form "<<lambda_H_max<<" to "<<
    //  lambda_H_min<< " of length "<< num_lambda_H<<std::endl;


    // @A: tau_max and tau_min are just auxiliary variables.

    double tau_H_max = log10(lambda_H_max);
    double tau_H_min = log10(lambda_H_min);

    arma::vec tau_H_vect(num_lambda_H);
    for (unsigned int i=0; i<num_lambda_H; i++)
    {
      tau_H_vect(i)=tau_H_max - i*(tau_H_max-tau_H_min)/(num_lambda_H -1);

      //  std::cout<<tau_vect(i)<<" "<<std::endl;
    }

    lambda_H = arma::exp10(tau_H_vect);

    // uncomment these lines to define an equispaced vector

    //     for (unsigned int i=0; i<num_lambda_H; i++)
    //     {
    //       lambda(i) = lambda_H_max - i*(lambda_H_max-lambda_H_min)/(num_lambda_H-1);
    //     }

    return lambda_H;
  }

  /*###########################################################################################
#
# @A Part5: Defining main Function, called "estimate", for definition_beta
#
###########################################################################################*/

  // MAIN FUNCTION TO COMPUTE THE ESTIMATION,
  // called from define_beta

  // @A: we don't have B and X here. They were defined as private variables. They are accessable when we are in the class.
  Rcpp::List estimate( arma::vec &tau, arma::vec &omega_vect,
                       // const int &N_iter,
                       Rcpp::Function computation_norm, Rcpp::Function Proc_Time,
                       double &thres, unsigned int &non_zeros_pred,
                       arma::vec &lambda_H, arma::vec &lambda_K,
                       bool verbose, Rcpp::Function sum_cpp)
  {

    /*###########################################################################################
#
# @A Part5-1: Introducing some new variables
#
###########################################################################################*/




    arma::mat B_subset(B.n_rows, B.n_cols-1);
    arma::mat X_subset(X.n_rows, X.n_cols-1);

    arma::vec X_j(num_data);
    arma::vec B_j(B.n_rows);

    arma::mat E_glob(E.n_rows, E.n_cols);
    arma::mat B_old(B.n_rows, B.n_cols);

    double lambda_H_iter;
    double lambda_K_iter;

    arma::vec vector_non_zeros_pred(num_pred);
    double num_zeros_pred;

    // @A: n_elem gives the number of elements in vector lambda
    int number_lambda_H_computed=lambda_H.n_elem;


    /*###########################################################################################
#
# @A Part5-2: Estimating Beta (iteration loops (3rd and 4th loops))
#
###########################################################################################*/

    // @A: In this part, there are 4 diffeent loops hirercial. lambda_K*lambda_H*N_iter*num_pred

    // @A: we have 3 different stopping criteria:
    // 1. N_iter
    // 2. \norm{B_old-B_new}_K < thres
    // 3. number of non zero predictors > non_zero_preds

    // @A: in this loop we don't care what the best lambda is. We just use diffenet lambdas from the vecotr lambda_H and lambda_K
    // we can find the best lambda_H and lambda_K based on BIC criteria or CV.


    arma::mat Time_mat(lambda_H.n_elem,lambda_K.n_elem);
    Time_mat.zeros();

    double remaining_time=0;
    bool BREAK_thres;

    arma::mat error_beta_matrix;

    double N_ite;

    for (unsigned int q=0; q<lambda_K.n_elem; q++) //loop on lambda_K
    {

      B.zeros();

      lambda_K_iter=lambda_K(q);

      for (unsigned int l=0; l<lambda_H.n_elem; l++) //loop on lambda_H
      {


        N_ite=0;

        BREAK_thres = TRUE;

        Time_mat(l,q)=Time(Proc_Time);

        lambda_H_iter=lambda_H(l);


        //@A: std::endl : means /n or break the line

        if(verbose) {Rcpp::Rcout<< std::endl<< "Lambda_H = "<< lambda_H_iter<< std::endl;}

        //      for (unsigned int k=0; k<N_iter; k++) // 3rd loop: B^{t} update for all N_iter
        do
        {

          // @A: in output, just print "*" between the results
          // number of * in output indicates the number of iterations.

          if(verbose) {Rcpp::Rcout << "*";}

          //  if (k==0 && l==0)
          // {
          //         E_glob=E;
          // }else{
          E_glob=Y-B*X.t(); //initialization with the previous estimation

          // }

          B_old=B;


          num_zeros_pred=0;
          vector_non_zeros_pred.zeros(num_pred,1);

          for (unsigned int j=0; j<num_pred; j++) // 4th loop: coordinate update
          {

            double omega=omega_vect(j);

            X_j=X.col(j);
            B_j= B_old.col(j);

            //if(k != 0)
            //    {

            // @A: B_j * X_j.t() is a k*N matrix. the first column is [B_1jX_j1,...,B_kjX_j1]' and etc. This is what we want for defining E.
            E= E_glob + B_j * X_j.t();

            // @A: B_tilde is the B_j^{v} in FLAME.
            arma::vec B_tilde = E * X_j / num_data;
            //  }else
            //    {
            //      B_tilde = B_j;
            //    }


            if (norm_H_vec(B_tilde) <= lambda_H_iter*omega)
            {
              B_j.zeros();
              num_zeros_pred ++;

            }else
            {

              // @A: "vector_non_zeros_pred" will become a p*1 vector which the first part of
              // them are nonzero (each element is the index of the predictor)
              // and rest of them are zero.
              vector_non_zeros_pred(j-num_zeros_pred)=j;

              double norm_B_j;

              // @A: computation_norm_cpp is used to calculate numerically the value of \norm{\beta_j}_H from the equation in AFSSEN
              norm_B_j = computation_norm_cpp(lambda_H_iter,lambda_K_iter, omega, B_tilde,
                                              tau, computation_norm); // numerical optimization to compute the norm_H

              // @A: This part shows why we put the projections of B and Y in the matrices.
              // Because K(B_j)=tau%B_j   (A%B=a_1*b_1+...+a_p*b_p)
              // The following come from the formula in AFSSEN by
              // B_j=((\norm{B_j}_H + \lambda_2 \omega_j)K + 2*\norm{B_j}_H*\lambda_1 L )^{-1} \norm{B_j}_H K(B_tilde)
              // which K(B_tilde)=\tau*B_tilde and if (K1+K2)(B)=cB then (K1+K2)^{-1}(B)=B/c
              B_j = (tau % B_tilde) * norm_B_j / (tau *(norm_B_j+lambda_H_iter*omega) + 2*lambda_K_iter*norm_B_j);

            }

            B.col(j) = B_j; // updating the B matrix

          } //end of the predictors loop

          // @A: end of the 4th loop

          error_beta_matrix = B_old-B;
          // check on the improvment of the estimation

          // Rcpp::Rcout << "number zeros pred = "<< num_zeros_pred << std::endl;

          // ?@A: error_beta_matrix may be better to calculate in norm H not norm K (norm K is slower).
          // if (norm_H_matrix(error_beta_matrix) < thres)
          // {
          // std::cout<<std::endl<<"Not significant change of Beta => Stop iteration at iter = "<<k+1<<std::endl;

          //@A: This loop is going to be run for N_iter times (N_iter is one of the stopping criteri).
          // But we have a stronger stopping criteria, a threshold on increment of estimated beta.
          // when it happens, we break the loop
          //    k=N_iter;

          //   if(verbose) {
          //   if(BREAK_thres){
          //   Rcpp::Rcout<< std::endl<< "BREAK based on Beta increament thresold" << std::endl;
          //   BREAK_thres=FALSE;
          //                   }
          //               }
          // }

          N_ite ++;

        } while(norm_H_matrix(error_beta_matrix) >= thres) ;


        if(verbose) {Rcpp::Rcout<< std::endl<< "number of Iteration = "<< N_ite << std::endl;}


        //end of the N_iter loop

        // @A: end of the 3rd loop

        /*###########################################################################################
#
# @A Part5-3: Estimating Beta ( lambda_K and lambda_H loops (1st and 2nd loops))
#
###########################################################################################*/

        vector_non_zeros_pred.resize(num_pred-num_zeros_pred);

        if (verbose)
        {
          Rcpp::Rcout<<"number of non zero predictors fitted with lambda_H = "<<lambda_H_iter<<" and lambda_K = "<<lambda_K_iter<<
            " is " <<num_pred-num_zeros_pred<<std::endl;
        }


        Time_mat(l,q)=Time(Proc_Time)-Time_mat(l,q);
        remaining_time=remaining_time+Time_mat(l,q);


        if( num_pred-num_zeros_pred >= non_zeros_pred )
        {

          // std::cout << "Number of non zeros predictors reached!";

          // @A: the +1 here is because R codifies c++ + 1.

          // A vector in R start from 1 bit in c++ start with 0

          number_lambda_H_computed=l+1;
          // std::cout<<"number lamba analyzed "<<number_lambda_computed<<std::endl;

          // @A: Another stopping criteria is a maximum for number of nonzero predictors, which is non_zeros_pred
          // if we reach to that number, we break the loop by following statement.
          l = lambda_H.n_elem ;
          if(q==0){
            if(verbose)
              std::cout <<std::endl<< "### It took = " << remaining_time << " seconds for a fixed lambda_K and " << number_lambda_H_computed<< " different lambda_H" <<std::endl;
          }


        }

        if((l==lambda_H.n_elem-1)&&(q==0)){

          if(verbose)
            std::cout <<std::endl<<std::endl<< "### It took = " << remaining_time << " seconds for a fixed lambda_K and " << number_lambda_H_computed<< " different lambda_H" <<std::endl;

        }


      }   //end of the lambda_H loop (l)
      // @A: end of the 2nd loop


    } //end of the lambda_K loop (q)

    // @A: end of the 1st loop

    /*###########################################################################################
#
# @A Part5-4: Return values ( and Finding the best lambda and estimate by BIC (it has not been activated in this code))
#
###########################################################################################*/

    //   std::cout<<"last value of lambda "<< lambda_iter<<" which identifies "<< num_zeros_pred<< " zeros predictors "<<std::endl;


    // @A: lambda_H_subset is a vector including the lambda_Hs were run in code
    //     (before breaking the loop by stopping criteria "non_zeros_pred")

    arma::vec lambda_H_subset=lambda_H.rows(0,number_lambda_H_computed-1);


    // @A: Actually the above loops will be break in a lambda and beta. We return them as final results.
    return(Rcpp::List::create(
        Rcpp::_["Beta"] = B,
        Rcpp::_["Pred"] = vector_non_zeros_pred+1,
        Rcpp::_["estimated_lambda_H"] = lambda_H_iter,
        Rcpp::_["Time_each_loop"] = Time_mat,
        Rcpp::_["Duration_Run_Code"] = sum_cpp(Time_mat),
        Rcpp::_["Lambda_H_vect"] = lambda_H_subset,
        Rcpp::_["lambda_H"] = lambda_H,
        Rcpp::_["lambda_K"] = lambda_K
    )
    );


  }

  /*###########################################################################################
#
# @A Part6: Cross Validation estimation
#
###########################################################################################*/

  // @A: Pay attention the X, Y and B are not input values of the following "estimate_CV" but we use them here.
  // It can be done because X, Y and B are the private values and can be accessable.


  Rcpp::List estimate_CV(const arma::mat &X_test, const arma::mat &Y_test,
                         arma::vec &tau, arma::vec &omega_vect,
                         Rcpp::Function computation_norm, Rcpp::Function Proc_Time,
                         double &thres, unsigned int &non_zeros_pred,
                         arma::vec &lambda_H, arma::vec &lambda_K,
                         double early_CV, double early_CV_thres, double max_ite, double target_inc,
                         bool verbose, Rcpp::Function sum_cpp)
  {


    arma::mat B_subset(B.n_rows, B.n_cols-1);
    arma::mat X_subset(X.n_rows, X.n_cols-1);

    arma::mat Beta_best;

    arma::vec X_j(num_data);
    arma::vec B_j(B.n_rows);

    arma::mat E_glob(E.n_rows, E.n_cols);
    arma::mat B_old(B.n_rows, B.n_cols);

    double lambda_H_iter;
    double lambda_K_iter;

    arma::vec vector_non_zeros_pred(num_pred);
    double num_zeros_pred;

    int number_lambda_H_computed=lambda_H.n_elem;
    int number_lambda_K_computed=lambda_K.n_elem;

    // @A: error_lambda_cv is gonna be a matrix including \norm{Y_estimated_for_test_data-Y_test}_K for each (lambda_H,lambda_K).
    arma::mat error_lambda_cv(number_lambda_H_computed,number_lambda_K_computed);
    arma::mat mat_target_inc(number_lambda_H_computed,number_lambda_K_computed);

    error_lambda_cv.zeros();
    mat_target_inc.zeros();

    double lambda_H_opt;
    double lambda_K_opt;

    double error_def;

    arma::vec vector_non_zeros_pred_opt;

    arma::vec vec_number_lambda_H(lambda_K.n_elem); // for each lambda_K, compute how many lambda_H we run before break


    arma::vec vec_break_in_number_pred(lambda_K.n_elem); // for each lambda_K, if kill switch is happened or not
    arma::vec vec_break_in_early_cv(lambda_K.n_elem); // for each lambda_K, if early CV is happened or not (base on earlyـCVـthreshold)

    arma::mat vec_norm_H_diff(lambda_H.n_elem,lambda_K.n_elem); // for each lambda_K, norm_H(B^{t+1} - B^{t})

    vec_norm_H_diff.zeros();

    arma::mat Time_mat_cv(lambda_H.n_elem,lambda_K.n_elem);
    Time_mat_cv.zeros();

    double remaining_time_cv=0;

    arma::mat error_beta_matrix;

    double N_iter_cv;

    arma::mat L_matrix(lambda_H.n_elem,max_ite);
    L_matrix.zeros();

    arma::mat Y_BX_ite;
    arma::vec B_ite;

    for (unsigned int q=0; q<lambda_K.n_elem; q++)
    {

      B.zeros();

      lambda_K_iter=lambda_K(q);

      vec_number_lambda_H(q)=lambda_H.n_elem;
      vec_break_in_number_pred(q)=0;
      vec_break_in_early_cv(q)=0;
      // @A: Pay attention: for first lambda_K we start with the B=0; and for each lambda_H we use a warm start
      // For next lambda_K we start with B=0


      for (unsigned int l=0; l<lambda_H.n_elem; l++)
      {

        Time_mat_cv(l,q)=Time(Proc_Time);

        lambda_H_iter=lambda_H(l);

        N_iter_cv=0;

        if (verbose) {Rcpp::Rcout<<std::endl<<"Lambda_H = "<< lambda_H_iter<< std::endl;}


        double aux_target=0;

        do{ // update for vector matrix beta^{t}
          if (verbose) {Rcpp::Rcout << "*";}

          B_old=B;

          num_zeros_pred=0;
          vector_non_zeros_pred.zeros(num_pred,1);

          for (unsigned int j=0; j<num_pred; j++)
          {

            double omega=omega_vect(j);
            E_glob = Y - B*X.t();

            X_j=X.col(j);
            B_j= B.col(j);

            E = E_glob + B_j * X_j.t();


            // E and E_glob are K*N matrices
            // X_j is a N*1 matrix
            // B_tilde is a K*1 matrix including the coefficient of B_tilde in basis

            arma::vec B_tilde = E * X_j / (num_data);


            if (norm_H_vec(B_tilde) <= lambda_H_iter*omega)
            {
              B_j.zeros();
              num_zeros_pred ++;

            }else
            {
              vector_non_zeros_pred(j-num_zeros_pred)=j;

              double norm_B_j;
              // in C++, A*B means matrix multiplication ... if A and B are m*N and N*s, the A*B will be m*s
              // in C++, A%B means elemntwise multiplication...if A and B are both N*1, the A%B will be N*1 by (a1*b1,...,b_N*b_N)


              norm_B_j = computation_norm_cpp(lambda_H_iter,lambda_K_iter, omega, B_tilde, tau, computation_norm);
              B_j = tau % B_tilde * norm_B_j / (tau*(norm_B_j+lambda_H_iter*omega) + lambda_K_iter*norm_B_j);

            }

            B.col(j) = B_j;


          } //end of the predictors loop
          // end of the 4th loop (coordinate update)


          //if(q==0){

          L_matrix(l,N_iter_cv)=target_function(B,tau,lambda_K_iter,lambda_H_iter,omega_vect);

          // if(verbose) {Rcpp::Rcout << std::endl<< "The value of L = "<<   L_matrix(l,N_iter_cv)  << std::endl;}

          //}


          if((target_inc==1) && (N_iter_cv > 0)){
            if(L_matrix(l,N_iter_cv) >= L_matrix(l,N_iter_cv-1)){
              aux_target=1;
            }
          }



          error_beta_matrix = B_old-B;

          N_iter_cv ++;

        } while( !(aux_target==1) && !(N_iter_cv == max_ite) && (norm_H_matrix(error_beta_matrix) > thres)); //end of the N_iter loop


        mat_target_inc(l,q)=aux_target;
        vec_norm_H_diff(l,q)=norm_H_matrix(error_beta_matrix);

        //           if(verbose) {Rcpp::Rcout << std::endl<< "number of Iteration = "<<   N_iter_cv  << std::endl;}

        // end of the 2nd loop (After each update of B)

        vector_non_zeros_pred.resize(num_pred-num_zeros_pred);

        //           if (verbose){
        //              Rcpp::Rcout<<"number of non zero predictors fitted with lambda_H = "<<lambda_H_iter<<" and lambda_K = "<<lambda_K_iter<<
        //                    " is " <<num_pred-num_zeros_pred<<std::endl;
        //          }


        arma::mat Y_pred_subset = B * X_test.t();
        arma::mat difference = Y_pred_subset-Y_test;

        // @A: We should calculate square_norm_H_matrix(difference). Y may not exist in K space.
        error_lambda_cv(l,q) = norm_H_matrix(difference);

        if((l==0) && (q==0)){
          error_def = error_lambda_cv(l,q)+1;
        }

        //           if(verbose) {Rcpp::Rcout<<"Cross validation is "<< error_lambda_cv(l,q) <<std::endl;}


        // For updating Beta_best

        //   std::cout<<error_lambda(l)<<error_lambda(l-1);
        if (error_lambda_cv(l,q) < error_def)
        {
          Beta_best=B;
          error_def = error_lambda_cv(l,q);
          lambda_H_opt=lambda_H(l);
          lambda_K_opt=lambda_K(q);
          vector_non_zeros_pred_opt=vector_non_zeros_pred;
        }



        //        if (verbose){
        //
        //             Rcpp::Rcout<<"error_def "<< error_def<<std::endl;
        //           Rcpp::Rcout<<"non zero predictors are "<< vector_non_zeros_pred + 1 ;
        //         Rcpp::Rcout<<"non zero predictors opt are "<< vector_non_zeros_pred_opt + 1 <<std::endl;
        //   }


        Time_mat_cv(l,q)=Time(Proc_Time)-Time_mat_cv(l,q);
        remaining_time_cv=remaining_time_cv+Time_mat_cv(l,q);



        //    if((l!=0) && verbose)  {Rcpp::Rcout<<"for l = "<< l << " and q = "<< q << " the |Diff(n,n-1)/cv(n-1)| will be " << fabs((error_lambda_cv(l-1,q)-error_lambda_cv(l,q))/error_lambda_cv(l-1,q)) <<std::endl;}

        // early_CV part
        // fabs returns absolute value for a flotting point

        if((l!=0) && (early_CV==1) && (fabs((error_lambda_cv(l-1,q)-error_lambda_cv(l,q))/error_lambda_cv(l-1,q)) < early_CV_thres) ){

          number_lambda_H_computed=l+1;
          l=lambda_H.n_elem;
          vec_number_lambda_H(q)=number_lambda_H_computed;
          vec_break_in_early_cv(q)=1;
          if(q==0){
            if(verbose)
              std::cout<<std::endl<<std::endl<< "### It took = " << remaining_time_cv << " seconds for a fixed lambda_K and " << number_lambda_H_computed<< " different lambda_H" <<std::endl;
          }
        }




        // kill_switch part
        if( num_pred-num_zeros_pred > non_zeros_pred )
        {
          // std::cout << "Number of non zeros predictors reached!";

          // @A: the +1 here is because R codifies c++ + 1.

          // A vector in R start from 1 bit in c++ start with 0

          number_lambda_H_computed=l+1;
          // std::cout<<"number lamba analyzed "<<number_lambda_computed<<std::endl;

          // @A: Another stopping criteria is a maximum for number of nonzero predictors, which is non_zeros_pred
          // if we reach to that number, we break the loop by following statement.
          l=lambda_H.n_elem ;

          vec_number_lambda_H(q)=number_lambda_H_computed;
          vec_break_in_number_pred(q)=1;


          if(q==0){

            if(verbose)
              std::cout<<std::endl<<std::endl<< "### It took = " << remaining_time_cv << " seconds for a fixed lambda_K and " << number_lambda_H_computed<< " different lambda_H" <<std::endl;

          }

        }


        if((l==lambda_H.n_elem-1)&&(q==0)){

          if(verbose)
            std::cout <<std::endl<<std::endl<< "### It took = " << remaining_time_cv << " seconds for a fixed lambda_K and " << number_lambda_H_computed<< " different lambda_H" <<std::endl;
        }



      } //end of the lambda_H loop

      // end of the 2nd loop


    } //end of the lambda_K loop

    // end of the 1st loop

    return Rcpp::List::create(
      Rcpp::_["Beta_opt"] = Beta_best,
      Rcpp::_["Pred_opt"] = vector_non_zeros_pred_opt+1,
      Rcpp::_["error"] = error_lambda_cv,
      Rcpp::_["Time"] = Time_mat_cv,
      Rcpp::_["Duration"] = sum_cpp(Time_mat_cv),
      Rcpp::_["lambda_H"] = lambda_H,
      Rcpp::_["lambda_K"] = lambda_K,
      Rcpp::_["lambda_H_opt"] = lambda_H_opt,
      Rcpp::_["lambda_k_opt"] = lambda_K_opt,
      Rcpp::_["vec_norm_H_diff_beta_inc"] = vec_norm_H_diff,
      Rcpp::_["vec_number_lambda_H_computed"] = vec_number_lambda_H,
      Rcpp::_["vec_break_in_number_pred"] = vec_break_in_number_pred,
      Rcpp::_["vec_break_in_early_cv"] = vec_break_in_early_cv,
      Rcpp::_["L_matrix"] = L_matrix,
      Rcpp::_["Mat_target_inc"] = mat_target_inc
    );


  }


};

// function to compute the estimation of Beta.
// called from R and explained in the file FLAME.R

//## @A:There are just two output main functions: Defining definition_beta(.) and Defining definition_beta_CV(.)

/*#####################
#
# @A Part: definition_beta(.) function
#
#####################*/

// @A: The main part of this function is the "estimate" function.
// The only difference between "definition_beta" and "estimate" is that "definition_beta"
// can compute automatically the vetor of lambda_H

// @A: Don't remember the f function in compputation_norm is gonna be come from R


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List definition_beta(const arma::mat &X, const arma::mat &Y,
                           arma::vec &tau, arma::vec &omega_vect,
                           Rcpp::Function computation_norm, Rcpp::Function Proc_Time,
                           double &thres, unsigned int &non_zeros_pred,
                           arma::vec &lambda_H_start, arma::vec &lambda_K_start,
                           double &ratio_lambda_H, int &num_lambda_H,
                           bool verbose,Rcpp::Function sum_cpp)
{
  linear_model LM(X,Y);

  arma::vec lambda_H_estimation;

  // @A: This part is for checking whether or not we manually insert the lambda vector.
  // If not, It automatically generate the log equspaced lambda_H with using the "definition_lambda_H" function.



  if (lambda_H_start.is_empty())
  {
    lambda_H_estimation=LM.definition_lambda_H(ratio_lambda_H, num_lambda_H, omega_vect, tau);
  }else{
    lambda_H_estimation=lambda_H_start;
  }
  //Rcpp::Rcout<<"lambda is not empty => ratio_lambda and num_lambda ignored!"<<std::endl;



  /*  if (Beta.n_cols==0)
  {
  LM.set_B_zeros();
  }
  */

  // @A: Define estimation as a Rcpp::list
  Rcpp::List estimation;
  // call the main function of the class linear_model
  estimation = LM.estimate(tau,  omega_vect,
                           computation_norm, Proc_Time,
                           thres, non_zeros_pred,
                           lambda_H_estimation, lambda_K_start,
                           verbose, sum_cpp);
  return(estimation);
}


// function to identify the best value of lambda with the cross validation method.
// called from R and described in FLAME.R
// N.B. the value of lambda_start MUST be a vector: the vector of all the
// values of lambda!!! Not empty!

/*#####################
#
# @A Part: Defining definition_beta_CV(.) function
#
#####################*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List definition_beta_CV(const arma::mat &X_train, const arma::mat &Y_train,
                              const arma::mat &X_test,  const arma::mat &Y_test,
                              arma::vec &tau, arma::vec &omega_vect,
                              Rcpp::Function computation_norm, Rcpp::Function Proc_Time,
                              double thres, unsigned int non_zeros_pred,
                              arma::vec &lambda_H_start, arma::vec &lambda_K_start,
                              double &ratio_lambda_H, int &num_lambda_H,
                              double early_CV, double early_CV_thres, double max_ite, double target_inc,
                              bool verbose, Rcpp::Function sum_cpp)
{


  // Here we initialized the X and Y with X_train and Y_train respectively. Actually Those
  // are used in estimate_CV part for computing the B.
  linear_model LM(X_train, Y_train);


  // @A: Unlike the "definition_beta" function, "definition_beta_CV" does not
  // generate the equispaced log lambda_H if we do not insert it.
  arma::vec lambda_H_estimation;

  lambda_H_estimation=lambda_H_start;
  if (lambda_H_start.is_empty())
  {
    lambda_H_estimation=LM.definition_lambda_H(ratio_lambda_H, num_lambda_H, omega_vect, tau);
  }

  // ?@A: Why not add a step for generating the equispaced in log for data
  // same as what we had in "definition_beta"?

  Rcpp::List est_CV;
  est_CV = LM.estimate_CV( X_test, Y_test,
                           tau, omega_vect,
                           computation_norm, Proc_Time,
                           thres, non_zeros_pred,
                           lambda_H_estimation, lambda_K_start,
                           early_CV, early_CV_thres, max_ite, target_inc,
                           verbose, sum_cpp);
  return(est_CV);
};


//End


