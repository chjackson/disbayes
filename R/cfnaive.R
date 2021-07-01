##' Simplistic estimator for case fatality given incidence, prevalence and mortality
##'
##' Estimator for case fatality given incidence, prevalence and mortality, under a simple discrete-time approximation. 
##'
##' @param inc Vector of estimates of the incidence rate by age
##'
##' @param prev Vector of estimates of the prevalence probability by age
##'
##' @param mort Vector of estimates of the mortality rate by age
##'
##' @details This estimator is based on a half-year discrete-time approximation and assumes that transitions only happen at half years. 
##' 
##' Let
##'
##' p12 = half-year transition probability from no disease to disease
##' 
##' p23 = half-year transition probability from disease to death
##'
##' then the probability of death within a year (two half-year cycles) averaging over people with and without the disease, 
##'
##' For a person who doesn't have the disease at the start, the probability that they die within a year is p12*p23.  For someone who has the disease at the start, they either die at the first cycle (probability p23) or survive the first cycle and die at the second (probability (1-p23)*p23)).   
##' 
##' Therefore the annual probability of death is computed by adding these three possibilities, weighted by the prevalence of disease at the start of the cycle. 
##'
##' mort = (1-prev)*p12*p23 + prev*(p23 + (1-p23)*p23)
##'
##' Therefore a compatible value for p23 can be obtained from
##' published prevalence, incidence and mortality as the solution for
##' x of the quadratic equation
##'
##' a x^2 + b*x + c = 0
##'
##' with a = -prev, b = (1-prev)*inc + 2*prev, c=-mort 
##' 
##' assuming the published incidence rate (inc) equals the probability
##' p12, and the case fatality rate (cf) equals the corresponding
##' probability p23, and the published mortality rate equals the
##' probability (mort) defined above. 
##'
##' @export
cfnaive <- function(inc, prev, mort){
    a <- prev
    b <- (1 - prev)*inc + 2*prev
    c <- - mort
    cf <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    cf
}
