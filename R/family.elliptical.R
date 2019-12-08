#' @title Family Objects for Elliptical Models
#' @method family elliptical
#' @name   family.elliptical
#' @aliases Cauchy
#' @aliases Cnormal
#' @aliases Gstudent
#' @aliases Glogis
#' @aliases Normal
#' @aliases Powerexp
#' @aliases Student
#' @aliases LogisI
#' @aliases LogisII
#' @description Family object provide an specify the details of the model used by functions such as \code{elliptical} and \code{gwer}. See the documentation for \code{elliptical} or \code{gwer} for more details.
#' @param object an object with the result of the fitted elliptical regression model.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return An object of class \dQuote{family} specifying a list with the follows elements:
#' \item{family}{character: the family name.}
#' \item{g0, g1, g2, g3, g4, g5}{derived fuctions associated with the distribution family defined.}
#' \item{df}{degree of freedom for t-student distribution.}
#' \item{s, r}{shape parameters for generalized t-student distribution.}
#' \item{alpha}{shape parameter for contaminated normal and generalized logistic distributions.}  
#' \item{mp}{shape parameter for generalized logistic distribution.}
#' \item{epsi,sigmap}{dispersion parameters for contaminated normal distribution.}
#' \item{k}{shape parameter for power exponential distribution.}
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{elliptical}}, \code{\link{gwer}}
#' @keywords Elliptical distributions
#' @keywords Elliptical models
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Normal()
#' ,data=luz)
#' family(elliptical.fitt)
#' @rdname family.elliptical
#' @export

family.elliptical<- function (object,...) 
{
  if (length(object$call$family) > 1) 
    eval(object$call$family)
  else eval(parse(text = paste(object$call$family, 
                               "()")))
}

#' @rdname family.elliptical
#' @method print family.elliptical
#' @noRd
#' @export
print.family.elliptical <- function (x, ...) 
{
  cat("\n", x$family, "family\n")
  cat("\n density  : ", as.character(as.list(x[["g0"]])), 
      "\n")
  cat("\n Wg: ", as.character(as.list(x[["g1"]])), "\n")
  cat("\n scale: ", as.character(as.list(x[["g2"]])), "\n")
  cat("\n scale dispersion: ", as.character(as.list(x[["g3"]])), 
      "\n")
  cat("\n scale variance: ", as.character(as.list(x[["g4"]])), 
      "\n")
  cat("\n Wg': ", as.character(as.list(x[["g5"]])), "\n")
  if (charmatch(x$family, "Gstudent", F)) 
    cat("\n r :", x$r, "\n", "\n s :", x$s, "\n")
  if (charmatch(x$family, "Glogis", F)) 
    cat("\n alpha :", x$alpha, "\n", "\n m :", x$mp, "\n")
  if (charmatch(x$family, "Cnormal", F)) 
    cat("\n epsilon :", x$epsi, "\n", "\n sigma :", x$sigmap, 
        "\n")
  if (charmatch(x$family, "Student", F)) 
    cat("\n df :", x$df, "\n")
  if (charmatch(x$family, "Powerexp", F)) 
    cat("\n k  :", x$k, "\n")
  return(invisible(0))
}



make.family.elliptical <- function (name, arg) 
{ 
  
  elliptical.deriv <- structure(
    .Data=list(
      g0=function(z,...) log(1/(sqrt(2*pi))*exp(-0.5*z^2)),
      g1=function(z,...) rep(-0.5,length(z)),
      g2=function(z,...) 1/4,
      g3=function(z,...) 3/4,
      g4=function(z,...) 1,
      g5=function(z,...) rep(0,length(z)),		   
      g0=function(z,...) log((1/pi)*(1+z^2)^(-1)),
      g1=function(z,...) -1/(1+z^2),
      g2=function(z,...) 1/8,
      g3=function(z,...) 3/8,
      g4=function(z,...) 1, ##non exist
      g5=function(z,...) 1/((1+z^2)^2),		   
      g0=function(z,df,...) log(((gamma(0.5*(df+1)))/((pi*df)^0.5*
                                                        gamma(0.5*df)))*(1+z^2/df)^(-0.5*(df+1))),
      g1=function(z,df,...) (df+1)/(-2*(df+z^2)),
      g2=function(z,df,...) (df+1)/(4*(df+3)),
      g3=function(z,df,...)(3*(df+1))/(4*(df+3)),
      g4=function(z,df,...) df/(df-2),
      g5=function(z,df,...) (df+1)/(2*(df+z^2)^2),
      g0=function(z,s,r,...) log((s^(r/2)*gamma(0.5+r/2))/(gamma(0.5)*
                                                             gamma(r/2))*(s+z^2)^(-0.5*(r+1))),
      g1=function(z,s,r,...) (r+1)/(-2*(s+z^2)),
      g2=function(z,s,r,...) (r*(r+1))/(4*s*(r+3)),
      g3=function(z,s,r,...) (3*(r+1))/(4*(r+3)),
      g4=function(z,s,r,...) s/(r-2),
      g5=function(z,s,r,...) (r+1)/(2*(s+z^2)^2),
      g0=function(z,...) log(1.4843300029*exp(-z^2)/(1+exp(-z^2))^2),
      g1=function(z,...) -tanh(z^2/2),
      g2=function(z,...) 1.477240176/4,
      g3=function(z,...) 4.013783934/4,
      g4=function(z,...) 0.79569,
      g5=function(z,...) -0.5+0.5*tanh(0.5*z^2)^2,
      g0=function(z,...) log(exp(z)/(1+exp(z))^2),
      g1=function(z,...) (exp(z)-1)/(-2*(z*(1+exp(z)))),
      g2=function(z,...) 1/12,
      g3=function(z,...) 2.42996/4,
      g4=function(z,...) pi^2/3,
      g5=function(z,...) (2*exp(z)*z-exp(2*z)+1)/(4*z^3*(1+exp(z)^2)),
      g0=function(z,alpha,mp,...) log((alpha*gamma(mp+mp))/(gamma(mp)*gamma(mp))*
                                        (exp(alpha*z)/(1+exp(alpha*z))^2)^mp),
      g1=function(z,alpha,mp,...) alpha*mp*(exp(alpha*z)-1)/(-2*(z*(1+exp(alpha*z)))), 
      g2=function(z,alpha,mp,...) (alpha^2*mp^2)/(4*(2*mp+1)),
      g3=function(z,alpha,mp,...) (2*mp)*(2+mp^2*trigamma(mp))/(4*(2*mp+1)),		   
      g4=function(z,alpha,mp,...) 2*trigamma(mp),
      g5=function(z,alpha,mp,...) alpha*mp*(2*alpha*exp(alpha*z)*z-exp(2*alpha*z)+1)/(4*z^3*(1+exp(alpha*z)^2)),
      g0=function(z,epsi,sigmap,...) log( (1-epsi)*1/(sqrt(2*pi))*exp(-0.5*z^2)+
                                            epsi*1/(sqrt(2*pi)*sigmap)*exp(-0.5*z^2/sigmap^2)),
      g1=function(z,epsi,sigmap,...)((1-epsi)*exp(-z^2/2)+
                                       (epsi*(sigmap^2)^(-1.5)*exp(-z^2/(2*sigmap^2))))/
        ((-2)*((1-epsi)*exp(-z^2/2)+
                 (epsi*(sigmap^2)^(-0.5)*exp(-z^2/(2*sigmap^2))))),
      g2=function(z,epsi,sigmap,...)
      {
        NULL
      },
      g3=function(z,epsi,sigmap,...) NULL, 
      g4=function(z,epsi,sigmap,...) 1+epsi*(sigmap^2-1),
      g5=function(z,epsi,sigmap,...)
      {
        NULL
      },
      g0=function(z,k,...) log(1/(gamma(1+((1+k)/2))*2^(1+(1+k)/2))*exp(-0.5*(abs(z)^(2/(1+k))))),		   
      g1=function(z,k,...) 1/(-2*(1+k)*(z^2)^(k/(1+k))),
      g2=function(z,k,...) (gamma((3-k)/2))/(4*(2^(k-1)*(1+k)^2*gamma((k+1)/2))),
      g3=function(z,k,...)(k+3)/(4*(k+1)) ,
      g4=function(z,k,...) 2^(1+k)*(gamma(1.5*(k+1))/(gamma((k+1)/2))),
      g5=function(z,k,...) k/(2*(z^2)^((2*k+1)/(1+k))*((1+k)^2))
    ),
    .Dim=c(6,9),
    .Dimnames=list(c("g0","g1","g2","g3","g4","g5"),
                   c("Normal","Cauchy","Student","Gstudent",
                     "LogisI","LogisII","Glogis",
                     "Cnormal","Powerexp")))
  
  
  
  
  if (is.character(name) && charmatch(name, dimnames(elliptical.deriv)[[2]], 
                                      F)) {
    g0 <- elliptical.deriv[["g0", name]]
    g1 <- elliptical.deriv[["g1", name]]
    g2 <- elliptical.deriv[["g2", name]]
    g3 <- elliptical.deriv[["g3", name]]
    g4 <- elliptical.deriv[["g4", name]]
    g5 <- elliptical.deriv[["g5", name]]
  }
  else {
    obj.deriv <- eval(parse(text = paste(name, ".deriv", 
                                         sep = "")))
    g0 <- obj.deriv[["g0", name]]
    g1 <- obj.deriv[["g1", name]]
    g2 <- obj.deriv[["g2", name]]
    g3 <- obj.deriv[["g3", name]]
    g4 <- obj.deriv[["g4", name]]
    g5 <- obj.deriv[["g5", name]]
  }
  family <- list(g0 = g0, g1 = g1, g2 = g2, g3 = g3, g4 = g4, 
                 g5 = g5, df = if (charmatch(name, "Student", F)) arg, 
                 s = if (charmatch(name, "Gstudent", F)) arg$s, r = if (charmatch(name, 
                                                                                  "Gstudent", F)) arg$r, alpha = if (charmatch(name, 
                                                                                                                               "Glogis", F)) arg$alpha, mp = if (charmatch(name, 
                                                                                                                                                                           "Glogis", F)) arg$m, epsi = if (charmatch(name, 
                                                                                                                                                                                                                     "Cnormal", F)) arg$epsi, sigmap = if (charmatch(name, 
                                                                                                                                                                                                                                                                     "Cnormal", F)) arg$sigmap, k = if (charmatch(name, 
                                                                                                                                                                                                                                                                                                                  "Powerexp", F)) arg)
  names(family) <- c("g0", "g1", "g2", "g3", "g4", "g5", "df", 
                     "s", "r", "alpha", "mp", "epsi", "sigmap", "k")
  structure(.Data = c(list(family = name), family), class = c("family.elliptical", 
                                                              "family"))
}

# Logistic Type I #
dlogisI <- function (x, mean = 0, sd = 1) 
{
  z <- (x - mean)/sd
  ((1.484300029 * exp(-(z^2)))/(sd * (1 + exp(-(z^2)))^2))
}


plogisI <- function (x, mean = 0, sd = 1) 
{
  z <- (x - mean)/sd
  (1.484300029/(sd * (1 + exp(-(z^2)))))
}


rlogisI <- function (n, mean = 0, sd = 1) 
{
  if (is.na(n)) 
    return(NA)
  u <- runif(n, 0, 1)
  mean+sd*sqrt(log(u/(1 - u)))
}


# Logistic Type II #
dlogisII <- function (x, mean = 0, sd = 1) 
{
  z <- (x - mean)/sd
  (exp(z)/(sd * (1 + exp(z))^2))
}


plogisII <- function (q) 
{
  punif((exp(q)/(1 + exp(q))), 0, 1)
}


rlogisII <- function (n) 
{
  if (is.na(n)) 
    return(NA)
  u <- runif(n, 0, 1)
  log(u/(1 - u))
}




# Power Exponential #
dpowerexp <- function (x, k = 0.5, mean = 0, sd = 1) 
{
  z <- (x - mean)/sd
  1/(sd * (gamma(1 + ((1 + k)/2)) * 2^(1 + (1 + k)/2))) * 
    exp(-0.5 * (abs(z)^(2/(1 + k))))
}


ppowerexp <- function (q, k = 0.5) 
{
  r <- 2/(1 + k)
  punif(q^r/2, -1, 1) * pgamma(q, (1 + 1/r), 1)
}


rpowerexp <- function (n, k = 0.5) 
{
  if (is.na(n)) 
    return(NA)
  u <- runif(n, -1, 1)
  r <- 2/(1 + k)
  ff <- rgamma(n, (1 + 1/r), 1)
  (2 * ff)^(1/r) * u
}



# Student Generalized #
dgstudent <- function (x, s = 1, r = 2, mean = 0, sd = 1) 
{
  z <- (x - mean)/sd
  dens <- (1/sd)(s^(r/2) * gamma(0.5 + r/2))/(gamma(0.5) * 
                                                gamma(r/2)) * (s + z^2)^(-0.5 * (r + 1))
}


pgstudent <- function (q, s = 1, r = 2) 
{
  (1 - pinvgamma(1/sqrt(q), s/2, r/2)) * pnorm(q, 0, 1)
}


rgstudent <- function (n, s = 1, r = 2) 
{
  if (is.na(n)) 
    return(NA)
  z <- rnorm(n, 0, 1)
  v <- rinvgamma(n, r/2, s/2)
  v^(-0.5) * z
}


# Inverse Gamma #
pinvgamma <- function (x, s = 1/2, r = 1) 
{
  if (is.na(x
  )) 
    return(NA)
  1 - pgamma((1/x), s, r)
}


rinvgamma <- function (n, s = 1/2, r = 1) 
{
  if (is.na(n)) 
    return(NA)
  1/(rgamma(n, s, r))
}



# Contaminated Normal #
rcnormal <- function(n, mean = 0, sd1 = 1, sd2, prob)
{
  sd <- sample(c(sd1,sd2), n, replace=T, prob=c(1-prob,prob))
  rnorm(n, mean, sd)
}





#' @rdname family.elliptical
#' @export

Normal <- function () 
{
  make.family.elliptical("Normal")
}



#' @rdname family.elliptical
#' @export

Cauchy <- function () 
{
  make.family.elliptical("Cauchy")
}



#' @rdname family.elliptical
#' @export

LogisI <- function () 
{
  make.family.elliptical("LogisI")
}



#' @rdname family.elliptical
#' @export

LogisII <- function () 
{
  make.family.elliptical("LogisII")
}



#' @rdname family.elliptical
#' @param df degrees of freedom.
#' @export

Student <- function (df = stop("no df argument")) 
{
  if (df < 0) 
    stop(paste("allowed values for degrees of freedom positive"))
  make.family.elliptical("Student", arg = df)
}



#' @rdname family.elliptical
#' @param k shape parameter.
#' @export


Powerexp <- function (k = stop("no k argument")) 
{
  if (abs(k) > 1) 
    stop(paste("k must be (-1,1)"))
  make.family.elliptical("Powerexp", arg = k)
}



#' @rdname family.elliptical
#' @param parma parameter vector (alpha, m).
#' @export

Glogis <- function (parma = stop("no alpha=alpha(m) or m argument")) 
{
  if ((parma[1] <= 0) || (parma[2] <= 0)) 
    stop(paste("alpha=alpha(m) and m must be positive"))
  make.family.elliptical("Glogis", arg = list(alpha = parma[1], 
                                              mp = parma[2]))
}

#' @rdname family.elliptical
#' @param parm parameter vector (s, r) for this distribuition.
#' @export

Gstudent <- function (parm = stop("no s or r argument")) 
{
  if ((parm[1] <= 0) || (parm[2] <= 0)) 
    stop(paste("s and r must be positive"))
  make.family.elliptical("Gstudent", arg = list(s = parm[1], 
                                                r = parm[2]))
}


#' @rdname family.elliptical
#' @param parmt parameters vector (epsi, sigma).
#' @export

Cnormal <- function (parmt = stop("no epsi or sigma argument")) 
{
  stop(paste("not implement yet"))
  if ((parmt[1] < 0) || (parmt[1] > 1) || (parmt[2] <= 0)) 
    stop(paste("0<=epsilon<=1 and sigma must be positive"))
  make.family.elliptical("Cnormal", arg = list(epsi = parmt[1], 
                                               sigmap = parmt[2]))
}



