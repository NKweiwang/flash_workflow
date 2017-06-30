# ATM function
#' title ash type model for l
#'
#' description use ash type model to maxmization
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'  }
#' @param Y is data matrix
#' @param Ef is mean for the factor
#' @param Ef2 is second moment for the factor
#' @param sigmae2 is variance structure for noise matrix
#' @param col_var is just for column variance, column is for column, row is for row
#' @param nonnegative is flag for whether output nonnegative value
#' @param output is format of output, "mean" is for mean value, "matrix" is for flash data matrix in ash
#'
#' @keywords internal
#'

ATM_r1 = function(Y, Ef, Ef2,
                  sigmae2, col_var = "row",
                  nonnegative = FALSE, output = "mean",
                  partype = "constant",ash_para = list()){
  # this part is preparing betahat and sebeta for ash
  if(is.matrix(sigmae2)){
    # print("sigmae2 is a matrix")
    # this is for all variance are different
    sum_Ef2 = (1/sigmae2) %*% Ef2
    sum_Ef2 = as.vector(sum_Ef2)
    sebeta = sqrt(1/(sum_Ef2))
    betahat = as.vector( (Y/sigmae2) %*% Ef ) / (sum_Ef2)
    betahat=as.vector(betahat)
  } else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
    # this is for the non-constant variance in column
    # print("sigmae2 is vector")
    if(col_var == "row"){
      sum_Ef2 = (1/sigmae2) * Ef2
      sum_Ef2 = sum(sum_Ef2)
      sebeta = sqrt(1/(sum_Ef2))
      betahat = as.vector( Y %*% (Ef/sigmae2) ) / (sum_Ef2)
      betahat=as.vector(betahat)
    }else {
      # for column
      # here I have alread change the dimension by taking transpose to Y
      sum_Ef2 = sum(Ef2)
      sebeta = sqrt( sigmae2/(sum_Ef2) )
      betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
      betahat=as.vector(betahat)
    }
  }else{
    # for constant case
    # print("sigmae2 is a constant")
    sum_Ef2 = sum(Ef2)
    sebeta = sqrt(sigmae2/(sum_Ef2))
    betahat = (t(Ef) %*% t(Y)) / (sum_Ef2)
    betahat=as.vector(betahat)
  }
  # set ash default setting up
  ash_default = list(betahat = betahat, sebeta = sebeta,
                     method = "fdr", mixcompdist = "normal")
  # ATM update
  # decide the sign for output
  if(nonnegative){
    ash_default$mixcompdist = "+uniform"
  }
  # decide the output to decide the convergence criterion
  # if(output == "matrix"){
    ash_para$outputlevel = 5
    ash_para = modifyList(ash_default,ash_para)
    ATM = do.call(ashr::ash,ash_para)
    # ATM = ash(betahat, sebeta, method="fdr", mixcompdist=mixdist,outputlevel=4)
    Ef = ATM$flash_data$postmean
#    SDf = ATM$result$PosteriorSD
    Ef2 = ATM$flash_data$postmean2
    if(output == "matrix"){
      mat_post = list(comp_postprob=ATM$flash_data$comp_postprob,
                    comp_postmean=ATM$flash_data$comp_postmean,
                    comp_postmean2=ATM$flash_data$comp_postmean2)
    } else {
      mat_post=NULL
    }
    fit_g = ATM$flash_data$fitted_g
    return(list(Ef = Ef,
                Ef2 = Ef2,
                mat = mat_post,
                g = fit_g))
  # } else {
  #   ash_para = modifyList(ash_default,ash_para)
  #   ATM = do.call(ashr::ash,ash_para)
  #   # ATM = ash(betahat, sebeta, method="fdr", mixcompdist=mixdist)
  #   Ef = ATM$flash_data$postmean
  # #  SDf = ATM$result$PosteriorSD
  #   Ef2 = ATM$flash_data$postmean2
  #   return(list(Ef = Ef, Ef2 = Ef2))
  # }

}

#' title prior and posterior part in objective function
#'
#' description prior and posterior part in objective function
#'
#' @return PrioPost the value for the proir and posterior in objectice function
#'  \itemize{
#'   \item{\code{PrioPost}} {the value for the proir and posterior in objectice function}
#'   \item{\code{penalty}} {penalty term of the g_l or g_f}
#'  }
#' @param mat matrix of flash.data in ash output which is for posterior
#' @param fit_g in fitted.g in ash output which is for prior
#' @keywords internal
#'
# two parts for the likelihood
Fval = function(mat,fit_g,ash_para){
  if(is.null(ash_para$method)){
    # in this case the method is default "fdr"
    # make sure the first row included in the index
    nonzeroindex = unique(c(1,which(fit_g$pi != 0)))
  }else{
    # there is no zero in the shrink method
    nonzeroindex = c(which(fit_g$pi != 0))
  }
  #prepare for the posterior part
  mat_postmean = mat$comp_postmean[nonzeroindex,]
  mat_postmean2 = mat$comp_postmean2[nonzeroindex,]
  mat_postvar = mat_postmean2 - mat_postmean^2
  mat_postprob = mat$comp_postprob[nonzeroindex,]
  # figure our the dimension
  K = dim(mat_postprob)[1]
  N = dim(mat_postprob)[2]
  if(is.vector(mat_postprob)){
    K = 1
    N = length(mat_postprob)
  }
  # prepare for the prior
  prior_pi = fit_g$pi[nonzeroindex]
  prior_var = (fit_g$sd[nonzeroindex])^2
  mat_priorvar = matrix(rep(prior_var,N),ncol = N)
  mat_priorprob = matrix(rep(prior_pi,N),ncol = N)
  # to get the value
  varodds = mat_priorvar / mat_postvar
  if(is.null(ash_para$method)){
    # only in the fdr case, which has zero
    varodds[1,] = 1 # deal with 0/0
  }
  probodds = mat_priorprob / mat_postprob
  ssrodds = mat_postmean2 / mat_priorvar
  if(is.null(ash_para$method)){
    # only in the fdr case, which has zero
    ssrodds[1,] = 1  # deal with 0/0
  }
  priorpost = mat_postprob * (log(probodds) - (1/2)*log(varodds) - (1/2)* (ssrodds -1))
  # in case of mat_postprob = 0
  priorpost[which(mat_postprob< 1e-100)] = 0
  PrioPost = sum(priorpost)

  #######
  # now we try to get the penalty term
  # I use the default of setting for lambda_k in https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/biostatistics/18/2/10.1093_biostatistics_kxw041/1/kxw041_Supp.pdf?Expires=1498001808&Signature=dtxDDFmYZ1u1jFoXfNg4mqGQnTPkQQMFZgt2-FE9r8GcgRl4DdEaP0oI6kzvUjZB5H3hLZmvAZfuFkt3GxojQJ1eX0y7I9kMJC75VTrxw1Ym~psIe8rdU5Xw4S4KZm79o1cOXBzY1iAORFkVM4z59uP1T1ltFJDu5fhlTYErU3flbiE1ivTm-BYN5P7w8dex-R5MWn98NpMVRUsoHRgbe7FhwavdnQ8QAD7~q1F4DzMyZk09TeNdEc7GNdBKbSOoS9mGGnvrqXvPz1nn~QzJ2U2CMCf04GfsgkclyUnH0aD6Vh1lSAfoEV47ZnT-9~942~f8rMs30H2Kqji7BjlctA__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q
  # this is a simple setting for fdr method we use the default lambda 10 for shrink we use nothing
  if(is.null(ash_para$method)){
    # in this case the method is default "fdr"
    #h_value = (10 - 1) * log(fit_g$pi[1])
    if(is.null(ash_para$nullweight)){
      h_value = (10 - 1) * log(fit_g$pi[1])
    }else{
      h_value = (ash_para$nullweight - 1) * log(fit_g$pi[1])
    }
  }else{
    # we assume that the only alternative method is shrink for now
    h_value = 0
  }
  return(list(PrioPost = PrioPost, penalty = h_value))
}

#' title conditional likelihood
#'
#' description conditional likelihood in objective function
#'
#' @return c_lik for the onditional likelihood in objectice function
#' @param N  dimension of residual matrix
#' @param P  dimension of residual matrix
#' @param sigmae2_v  residual matrix
#' @param sigmae2  variance structure of error
#' @keywords internal
#'
# this version we need to know the truth of sigmae2, we can use sigmae2_v as the truth
C_likelihood = function(N,P,sigmae2_v,sigmae2){
  if(is.matrix(sigmae2)){
    # print("clik as matrix of sigmae2")
    c_lik = -(1/2) * sum( log(2*pi*sigmae2) + (sigmae2_v)/(sigmae2) )
  }else if(is.vector(sigmae2) & length(sigmae2)==P){
    # print("clik as vector sigmae2")
    # change the format to fit the conditional likelihood
    sigmae2_v = colMeans(sigmae2_v)
    c_lik = -(N/2) * sum( log(2*pi*sigmae2) + (sigmae2_v)/(sigmae2) )
  } else {
    # print("clik as sigmae2 is constant")
    # change the format to fit the conditional likelihood and accelerate the computation.
    sigmae2_v = mean(sigmae2_v)
    c_lik = -(N*P)/2 * ( log(2*pi*sigmae2) + (sigmae2_v)/(sigmae2) )
    # here I want to use fully variantional inference Elogsigmae2 rather logEsigmae2
    # c_lik = -(N*P)/2 * ( log(2*pi*sigmae2_true) + (sigmae2_v)/(sigmae2_true) + log(N*P/2) - digamma(N*P/2) )
  }
  return(list(c_lik = c_lik))
}

#' title objective function in VEM
#'
#' description  objective function
#'
#' @return obj_val value of the objectice function
#' @param N  dimension of residual matrix
#' @param P  dimension of residual matrix
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure (we use the estimated one to replace that if the truth is unknown)
#' @param par_l ash output for l
#' @param par_f ash output for f
#' @param objtype  objective function type,
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @keywords internal
#'

# objective function
obj = function(N,P,sigmae2_v,sigmae2,par_f,par_l,objtype = "margin_lik",ash_para){
  if(is.list(sigmae2)){
    # print("obj using kronecker product")
    sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
  }
  if(objtype=="lowerbound_lik"){
    priopost_f = Fval(par_f$mat, par_f$g,ash_para)$PrioPost
    priopost_l = Fval(par_l$mat, par_l$g,ash_para)$PrioPost
    penalty_f = Fval(par_f$mat, par_f$g,ash_para)$penalty
    penalty_l = Fval(par_l$mat, par_l$g,ash_para)$penalty
    c_lik = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
    obj_val = c_lik + priopost_l + priopost_f + penalty_l + penalty_f
  } else {
    obj_val = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
  }
  return(obj_val)
}

#' title rescale the lambda_l and lambda_f for identifiablity
#'
#' description rescale the lambda_l and lambda_f for identifiablity
#'
#' @return a list of sig2_l and sig2_f for the rescaled variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @keywords internal
#'
rescale_sigmae2 = function(sig2_l,sig2_f){
  norm_l = sqrt(sum(sig2_l^2))
  norm_f = sqrt(sum(sig2_f^2))
  norm_total = norm_l * norm_f
  sig2_l = sig2_l / norm_l
  sig2_f = sig2_f / norm_f
  sig2_l = sig2_l * sqrt(norm_total)
  sig2_f = sig2_f * sqrt(norm_total)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}

#' title initial value for Bayes variance structure estimation for kronecker productor
#'
#' description initial value for Bayes variance structure estimation for kronecker productor
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @keywords internal
#'
inital_Bayes_var = function(sigmae2_v){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  # estimate the initial value of lambda_l and lambda_f
  sig2_l_pre = rowMeans(sigmae2_v)
  sig2_f_pre = colMeans( sigmae2_v / matrix(rep(sig2_l_pre,P), ncol = P) )
  sig2_pre_list = rescale_sigmae2(sig2_l_pre,sig2_f_pre)
  sig2_l_pre = sig2_pre_list$sig2_l
  sig2_f_pre = sig2_pre_list$sig2_f
  # start the iteration
  maxiter = 100
  inital_tol = 1e-3
  tau = 0
  epsilon = 1
  while(epsilon > inital_tol & tau <= maxiter){
    tau = tau + 1
    sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) )
    sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) )
    sig2_list = rescale_sigmae2(sig2_l,sig2_f)
    sig2_l = sig2_list$sig2_l
    sig2_f = sig2_list$sig2_f
    epsilon = sqrt(mean((sig2_f - sig2_f_pre)^2 )) + sqrt(mean((sig2_l - sig2_l_pre)^2))
    sig2_l_pre = sig2_l
    sig2_f_pre = sig2_f
  }
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}

#' title Bayes variance structure estimation for kronecker productor
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2 variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'
Bayes_var = function(sigmae2_v,sigmae2 = NA){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  if( is.na(sigmae2) || !is.list(sigmae2) ){
    # we don't know the truth this is in the first iteration
    sigmae2 = inital_Bayes_var(sigmae2_v)
  }
  sig2_l_pre = sigmae2$sig2_l
  sig2_f_pre = sigmae2$sig2_f
  # this has already beed rescaled
  # here we use alpha_l = alpha_f = beta_l = beta_f = 0
  sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) )
  sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) )
  #rescaled the variance
  sig2_list = rescale_sigmae2(sig2_l,sig2_f)
  sig2_l = sig2_list$sig2_l
  sig2_f = sig2_list$sig2_f
  # sig2_out = matrix(rep(sig2_l,P),ncol = P) * matrix(rep(sig2_f,each = N),ncol = P)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}


#' title noisy variance structure estimation
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'

noisy_var = function(sigmae2_v,sigmae2_true){
  # a reaonable upper bound which control optimal algorithm
  upper_range = abs(mean(sigmae2_v - sigmae2_true)) + sd(sigmae2_v - sigmae2_true)
  # value of likelihood
  # f_lik <- function(x,sigmae2_v_l = sigmae2_v,sigmae2_true_l = sigmae2_true){
  #  -(1/2)*sum( log(2*pi*(x+sigmae2_true_l)) + (1/(x+sigmae2_true_l))*sigmae2_v_l )
  # }
  # I need the negative of the f_lik since optim find the minimum rather than the maximum
  f_lik <- function(x,sigmae2_v_l = sigmae2_v,sigmae2_true_l = sigmae2_true){
    (1/2)*mean( log(2*pi*(x+sigmae2_true_l)) + (1/(x+sigmae2_true_l))*sigmae2_v_l )
  }
  AA <- optim(mean(sigmae2_v - sigmae2_true), f_lik,lower = 0, upper = upper_range,  method = "Brent")
  e_sig = AA$par
  return(e_sig + sigmae2_true)
}


#' title noisy variance structure estimation with noisy variance on each column
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'

noisy_var_column = function(sigmae2_v,sigmae2_true){
  # in this case, we need to think about the for each column
  # we can use noisy_var = function(sigmae2_v,sigmae2_true) with the input of each column
  P = dim(sigmae2_true)[2]
  sig2_out = sapply(seq(1,P),function(x){noisy_var(sigmae2_v[,x],sigmae2_true[,x])})
  return(sig2_out)
}


#' title module for estiamtion of the variance structure
#'
#' description estiamtion of the variance structure
#'
#' @return sigmae2 estimated variance structure
#' @param partype parameter type for the variance,
#' "constant" for constant variance,
#' "var_col" for nonconstant variance for column,
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "noisy" is noisy version s_ij + sigmae2
#' "noisy_col" noisy version for column s_ij + sigmae2_j
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure
#' @param sigmae2_pre the previouse sigmae2 which is only used in Bayes_var
#' @keywords internal
#'
# sigma estimation function
sigma_est = function(sigmae2_v,sigmae2_true,sigmae2_pre = NA, partype = "constant"){
  if(partype == "var_col"){
    # print("use var col method")
    sigmae2 = colMeans(sigmae2_v)
  } else if(partype == "noisy"){
    # print("use the noisy version")
    sigmae2 = noisy_var(sigmae2_v,sigmae2_true)
  } else if(partype == "noisy_col"){
    sigmae2 = noisy_var_column(sigmae2_v,sigmae2_true)
  }else if (partype == "Bayes_var"){
    # here sigmae2_true is a list
    # print("use the kronecker product")
    sigmae2 = Bayes_var(sigmae2_v,sigmae2_pre)
  }else if (partype == "known"){
    # print("use the known sigmae")
    sigmae2 = sigmae2_true
  } else {
    # print("use the mean value for sigmae")
    # this is for constant case
    sigmae2 = mean(sigmae2_v)
  }
  return(sigmae2)
}

#' title empirical sample variance matrix estimation
#'
#' description this function is to calculate the empirical sample variance matrix sigmae2_v
#'
#' @return sigmae2_v  empirical sample variance matrix
#' @param Y which is the residual
#' @param fl_list this is a list for the El Ef El2 and Ef2 from other factors
#' @param El mean estimation for the current loading
#' @param Ef mean estimation for the current factor
#' @param El2 second momnet estimation for the current factor
#' @param Ef2 second moment estimation for the current loading
#' @keywords internal
#'
sigmae2_v_est = function(Y,El,Ef,El2,Ef2,fl_list=list()){
  # in this version of function we haven't included the missing value case in
  if(length(fl_list)==0){
    # this is for the rank one case since there are no other factors
    # here Y is just the original data
    sigmae2_v =  Y^2 - 2*Y*(El %*% t(Ef)) + (El2 %*% t(Ef2))
  }else{
    # this is the for the rank more than one we have other factors
    # now the Y is residual matrix rather the original matrix
    # to get the original data we use Yhat
    Yhat = Y + fl_list$El %*% t(fl_list$Ef)
    # the residual matrix should be
    # this is for E(l_1f_1+l_2f_2+...+l_kf_k)^2
    fl_norm = (El%*%t(Ef) + fl_list$El%*%t(fl_list$Ef))^2 -
      (El^2 %*% t(Ef^2) + (fl_list$El)^2 %*% t((fl_list$Ef)^2)) +
      (El2 %*% t(Ef2) + fl_list$El2 %*% t(fl_list$Ef2))
    sigmae2_v = Yhat^2 - 2* Yhat * (El%*%t(Ef) + fl_list$El%*%t(fl_list$Ef)) + fl_norm
  }
  return(sigmae2_v)
}


#' title one step update in flash iteration using ash
#'
#' description one step update in flash iteration using ash
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{Ef2}} {is a N vector for second moment of factors}
#'   \item{\code{sigmae2_v}}{is a N by P matrix for residual square}
#'   \item{\code{sigmae2_true}}{is a N by P matrix for estimated value for the variance structure}
#'   \item{\code{obj_val}}{the value of objectice function}
#'  }
#' @param Y the data matrix
#' @param N dimension of Y
#' @param P dimension of Y
#' @param El mean for the loadings
#' @param El2 second moment for the loadings
#' @param Ef mean for the factors
#' @param Ef2 second moment for the factors
#' @param sigmae2_v residual square
#' @param sigmae2_true the (true) known variance structure
#' Here, sigmae2 is the estimated variance structure in each step
#' sigmae2_true is the truth we know, some times sigmae2 is noisy version of sigmae2_true
#' @param sigmae2 the estimation of the variance structure
#' @param nonnegative if the facotor and loading are nonnegative or not.
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param partype parameter type for the variance,
#' "constant" for constant variance,
#' "var_col" for nonconstant variance for column,
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param objtype  objective function type,
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint
#' @keywords internal
#'
# one step update function
one_step_update = function(Y, El, El2, Ef, Ef2,
                           N, P,
                           sigmae2_v, sigmae2_true,
                           sigmae2,
                           nonnegative = FALSE,
                           partype = "constant",
                           objtype = "margin_lik",
                           fix_factor = FALSE,
                           ash_para = list(),
                           fl_list=list()){
  # if fix_factor is True, please choose objtype = "margin_lik"
  output = ifelse(objtype == "lowerbound_lik", "matrix", "mean")
  # deal with the missing value
  na_index_Y = is.na(Y)
  is_missing = any(na_index_Y)  # keep the missing result
  if(is_missing){
    # print("missing value use the EY as missing Y")
    Y[na_index_Y] = (El %*% t(Ef))[na_index_Y]
  }
  if(fix_factor){
    # actually we need do nothing here
    output = "mean"
    objtype = "margin_lik"
  }else{
    # estimate the variance structure
    # sigma_est = function(sigmae2_v,sigmae2_true,sigmae2_pre = NA, partype = "constant")
    sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
    # for the list case in kronecker product
    if(is.list(sigmae2)){
      # print("sigmae is kronecker product in estimating sigmae2 in one step update")
      sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
    }
    # Y = lf^T + E and ATM is for l given f, for f given l we need y^T = fl^T + E^T
    # sigmae2_input = ifelse(is.matrix(sigmae2),t(sigmae2),sigmae2) is wrong here
    if(is.matrix(sigmae2)){
      # print("sigmae is matrix in estimating sigmae2 in one step update")
      sigmae2_input = t(sigmae2)
    }else{
      sigmae2_input = sigmae2
    }
    par_f = ATM_r1(t(Y), El, El2, sigmae2_input,
                   col_var = "column", nonnegative,
                   output,partype,ash_para)
    Ef = par_f$Ef
    Ef2 = par_f$Ef2
    # if the Ef is zeros ,just return zeros
    if(sum(Ef^2) <= 1e-12){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      sigmae2_v =  sigmae2_v_est(Y,El,Ef,El^2,Ef^2,fl_list)
      sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
      obj_val = obj(N, P, sigmae2_v, sigmae2, par_f=NA, par_l=NA, objtype="margin_lik",ash_para)
      return(list(El = El, El2 = El^2,
                  Ef = Ef, Ef2 = Ef^2,
                  sigmae2_v = sigmae2_v,
                  sigmae2 = sigmae2,
                  obj_val = obj_val))
    }
    # update the Y and the Y^2 in this if else block
    if(is_missing){
      # the missing varianece is just sigmae2_v + sigmae2
      if(is.matrix(sigmae2)){
        sigmae2_impute = sigmae2[na_index_Y]
      }else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
        sigmae2_impute = (matrix(rep(sigmae2,each = N),ncol = P))[na_index_Y]
      }else{
        sigmae2_impute = sigmae2
      }
      sigmae2_v = sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
      sigmae2_v[na_index_Y] = sigmae2_v[na_index_Y] + sigmae2_impute
    }else{
      sigmae2_v = sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
    }
  }

  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
  if(is.list(sigmae2)){
    # kronecker product
    sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
  }
  par_l = ATM_r1(Y, Ef, Ef2,
                 sigmae2, col_var = "row",
                 nonnegative, output,
                 partype, ash_para)
  El = par_l$Ef
  El2 = par_l$Ef2
  # if El is zeros just return
  if(sum(El^2) <= 1e-12){
    El = rep(0,length(El))
    Ef = rep(0,length(Ef))
    sigmae2_v = sigmae2_v_est(Y,El,Ef,El^2,Ef^2,fl_list)
    sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
    obj_val = obj(N, P, sigmae2_v, sigmae2, par_f=NA, par_l=NA, objtype="margin_lik",ash_para)
    return(list(El = El, El2 = El^2,
                Ef = Ef, Ef2 = Ef^2,
                sigmae2_v = sigmae2_v,
                sigmae2 = sigmae2,
                obj_val = obj_val))
  }
  if(is_missing){
    # the missing varianece is just sigmae2_v + sigmae2
    if(is.matrix(sigmae2)){
      sigmae2_impute = sigmae2[na_index_Y]
    }else if(is.vector(sigmae2) & length(sigmae2) == length(Ef)){
      sigmae2_impute = (matrix(rep(sigmae2,each = N),ncol = P))[na_index_Y]
    }else{
      sigmae2_impute = sigmae2
    }
    sigmae2_v = sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
    sigmae2_v[na_index_Y] = sigmae2_v[na_index_Y] + sigmae2_impute
  }else{
    sigmae2_v =  sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
  }
  #use the estiamtion as the truth
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
  obj_val = obj(N, P, sigmae2_v, sigmae2, par_f, par_l, objtype,ash_para)

  return(list(El = El, El2 = El2,
              Ef = Ef, Ef2 = Ef2,
              sigmae2_v = sigmae2_v,
              sigmae2 = sigmae2,
              obj_val = obj_val))
}

#' inital value for flash
#'
#' description inital value for flash
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{El2}} {is a N vector for second moment of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{Ef2}} {is a N vector for second moment of factors}
#'   \item{\code{sigmae2_v}}{is a N by P matrix for residual square}
#'   \item{\code{sigmae2_true}}{is a N by P matrix for estimated value for the variance structure}
#'  }
#' @param Y the data matrix
#' @param nonnegative if the facotor and loading are nonnegative or not.
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint
#' @param factor_value is the factor value if the factor is fixed
#' @param fl_list this is a list for the El, El2, Ef and Ef2 for other factor loadings from the model.
#'                in the flash for backfitting, we need those to get the residuals and estimation for variance
#' @param initial_list_r1 this is a list for the inital values of current El Ef El2 and Ef2.
#' @param fix_initial this is a indicator for the initial value is fix or not.
#' @keywords internal
#'
initial_value = function(Y, nonnegative = FALSE,
                         factor_value = NA,fix_factor = FALSE,
                         initial_list_r1 = list(), fix_initial = FALSE,
                         fl_list=list()){
  # use the total mean as the estimated missing value
  na_index_Y = is.na(Y)
  is_missing = any(na_index_Y)
  if(is_missing){
    # this is the initialization for the missing value
    Y[na_index_Y] = mean(Y, na.rm = TRUE)
    # Y[na_index_Y] = 0
  }
  # the flash with plug-in missing value
  if(fix_initial){
    El = initial_list_r1$El
    Ef = initial_list_r1$Ef
    El2 = initial_list_r1$El2
    Ef2 = initial_list_r1$Ef2
  }else if(fix_factor){
    Ef = factor_value
    Ef2 = Ef^2
    El = as.vector( (Y %*% Ef) / (sum(Ef2)) )
    El2 = El^2
  }else{
    # El = svd(Y)$u[,1]
    El = as.vector(irlba::irlba(Y,nv = 0,nu = 1)$u)
    # the nonnegative value need positive inital value
    if(nonnegative){
      # El = abs(El)
      El = El * (El > 0)
    }
    El2 = El^2
    Ef = as.vector(t(El)%*%Y)
    Ef2 = Ef^2
  }
  # residual matrix initialization
  # here we don't have any information about the sigmae2_true
  if(is_missing){
    # here we don't have any information about the variance and the Y as well
    # since the previoust the Y is impute by the mean of observed Ys so I would use the mean of Y here as well.
    # this initialization still need to be considered more.
    sigmae2_v =  sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
  }else{
    sigmae2_v =  sigmae2_v_est(Y,El,Ef,El2,Ef2,fl_list)
  }
  return(list(El = El, El2 = El^2,
              Ef = Ef, Ef2 = Ef^2,
              sigmae2_v = sigmae2_v))
}


#' objective function value at zero
#'
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{obj_val}} {objective function value at zero}
#'  }
#' @param Y the data matrix
#' @param sigmae2_true true value for the variance structure
#' @param partype parameter type for the variance,
#' "constant" for constant variance,
#' "var_col" for nonconstant variance for column,
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param fl_list is a list containing all the informations from other factors ans loadings
#' l is loadings, f is factors, l2 and f2 are corresponding second moments.
#' priorpost_vec is expectation of log piropost ratio
#' clik is conditional likelihood (marginal likelihood)
#'
#' @details objective function value at zero, use to check the rank one estimation and rank zero estimation results.
#'

obj_zero = function(Y,N,P,sigmae2_true,fl_list,partype){
  na_index_Y = is.na(Y)
  is_missing = any(na_index_Y)
  if(is_missing){
    # this is the initialization for the missing value
    Y[na_index_Y] = mean(Y, na.rm = TRUE)
  }
  El = rep(0,N)
  Ef = rep(0,P)
  sigmae2_v = sigmae2_v_est(Y,El,Ef,El^2,Ef^2,fl_list)
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true, NA ,partype)
  obj_val = obj(N, P, sigmae2_v, sigmae2, par_f=NA, par_l=NA, objtype="margin_lik",ash_para = list(method = NULL))
  return(list(obj0 = obj_val,
              sigmae2_rank0_est = sigmae2))
}


#' FLASH
#'
#' factor loading adaptive shrinkage rank one version
#' @return list of factor, loading and variance of noise matrix
#'  \itemize{
#'   \item{\code{El}} {is a N vector for mean of loadings}
#'   \item{\code{Ef}} {is a N vector for mean of factors}
#'   \item{\code{sigmae2}}{is a N by P matrix for estimated value for the variance structure}
#'  }
#' @param Y the data matrix
#' @param tol is for the tolerence for convergence in iterations and ash
#' @param maciter_r1 is maximum of the iteration times for rank one case
#' @param sigmae2_true true value for the variance structure
#' @param nonnegative if the facotor and loading are nonnegative or not.
#' TRUE for nonnegative
#' FALSE for no constraint
#' @param partype parameter type for the variance,
#' "constant" for constant variance,
#' "var_col" for nonconstant variance for column,
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "loganova" is anova estiamtion for the log residual square
#' @param objtype  objective function type,
#' "margin_lik" for conditional likelihood,
#' "lowerbound_lik" for full objective function
#' @param fix_factor whether the factor is fixed or not
#' TRUE for fix_factor
#' FALSE for non-constraint
#' @param factor_value is the factor value if the factor is fixed
#' @param ash_para is the parameters list for ash
#' @param fl_list is a list containing all the informations from other factors ans loadings
#' l is loadings, f is factors, l2 and f2 are corresponding second moments.
#' priorpost_vec is expectation of log piropost ratio
#' clik is conditional likelihood (marginal likelihood)
#'
#' @details flash privide rank one matrix decomposition with variational EM algorithm.
#'
#' @export flash
#'
#' @importFrom ashr ash
#'
flash = function(Y, tol=1e-5, maxiter_r1 = 500,
                 partype = c("constant","known","Bayes_var","var_col","noisy","noisy_col"),
                 sigmae2_true = NA,
                 factor_value = NA,fix_factor = FALSE,
                 initial_list_r1 = list(), fix_initial = FALSE,
                 nonnegative = FALSE,
                 objtype = c("margin_lik","lowerbound_lik"),
                 ash_para = list(),
                 fl_list=list()){
  # match the parameters
  partype = match.arg(partype, c("constant","known","Bayes_var","var_col","noisy","noisy_col"))
  objtype = match.arg(objtype, c("margin_lik","lowerbound_lik"))

  # check the input
  if( any(!is.na(sigmae2_true)) & !(partype %in% c("known","noisy","noisy_col")) ){
    stop("You should choose partype = 'known' or 'noisy' ,'noisy_col',if you want to input the value of sigmae2_true")
  }

  N = dim(Y)[1]
  P = dim(Y)[2]

  # check the obj value at zero
  obj0_list = obj_zero(Y,N,P,sigmae2_true,fl_list,partype)

  # to get the inital values
  g_initial = initial_value(Y, nonnegative,factor_value,fix_factor,initial_list_r1, fix_initial,fl_list)
  El = g_initial$El
  Ef = g_initial$Ef
  El2 = g_initial$El2
  Ef2 = g_initial$Ef2
  sigmae2_v = g_initial$sigmae2_v

  # start iteration
  # in the first iteration, there is no value foe estimated sigmae2
  g_update = one_step_update(Y, El, El2, Ef, Ef2,
                             N, P,
                             sigmae2_v, sigmae2_true,
                             sigmae2 = NA,
                             nonnegative ,
                             partype ,
                             objtype ,
                             fix_factor,
                             ash_para,
                             fl_list)
  # parameters updates
  El = g_update$El
  El2 = g_update$El2
  Ef = g_update$Ef
  Ef2 = g_update$Ef2
  sigmae2_v = g_update$sigmae2_v
  sigmae2 = g_update$sigmae2
  obj_val = g_update$obj_val

  # track the objective value
  obj_val_track = c(obj_val)

  ##########
  print(sigmae2)

  # we should also return when the first run get all zeros
  if(sum(El^2)==0 || sum(Ef^2)==0){
    sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
    if(is.list(sigmae2)){
      # print("here using kronecker product")
      sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
    }
    # add one more output for greedy algorithm which not useful here
    c_lik_val = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
    # the above value is not useful, but is helpful to get the postprior value
    # since obj_val = c_lik_value + priorpost_l + priorpost_f
    return(list(l = El, f = Ef, l2 = El2, f2 = Ef2,
                sigmae2 = sigmae2,
                obj_val = obj_val,
                c_lik_val = c_lik_val))
  }

  epsilon = 1
  tau = 1
  while(epsilon >= tol & tau < maxiter_r1){
    tau = tau + 1
    pre_obj = obj_val

    g_update = one_step_update(Y, El, El2, Ef, Ef2,
                               N, P,
                               sigmae2_v, sigmae2_true,
                               sigmae2,
                               nonnegative ,
                               partype ,
                               objtype ,
                               fix_factor,
                               ash_para,
                               fl_list)
    # parameters updates
    El = g_update$El
    El2 = g_update$El2
    Ef = g_update$Ef
    Ef2 = g_update$Ef2
    sigmae2_v = g_update$sigmae2_v
    sigmae2 = g_update$sigmae2
    obj_val = g_update$obj_val


    ##########
    print(sigmae2)

    if(sum(El^2)==0 || sum(Ef^2)==0){
      El = rep(0,length(El))
      Ef = rep(0,length(Ef))
      break
    }
    epsilon = abs(pre_obj - obj_val)
    obj_val_track = c(obj_val_track,obj_val)
  }
  sigmae2 = sigma_est(sigmae2_v,sigmae2_true,sigmae2 ,partype)
  if(is.list(sigmae2)){
    # print("here, using kronecker product")
    sigmae2 = sigmae2$sig2_l %*% t(sigmae2$sig2_f)
  }
  # add one more output for greedy algorithm which not useful here
  c_lik_val = C_likelihood(N,P,sigmae2_v,sigmae2)$c_lik
  # the above value is not useful, but is helpful to get the postprior value
  # since obj_val = c_lik_value + priorpost_l + priorpost_f

  # add a check step
  # if(obj_val < obj0_list$obj0){
  #   El = rep(0,N)
  #   Ef = rep(0,P)
  #   return(list(l = El, f = Ef, l2 = El^2, f2 = Ef^2,
  #               sigmae2 = obj0_list$sigmae2_rank0_est,
  #               obj_val = obj0_list$obj0,
  #               c_lik_val = obj0_list$obj0))
  # }

  return(list(l = El, f = Ef, l2 = El2, f2 = Ef2,
              sigmae2 = sigmae2,
              obj_val = obj_val,
              c_lik_val = c_lik_val,
              obj_val_track = obj_val_track))
}
