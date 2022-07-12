state_probs_norem <- function(i, f)  {
    l <- i + f
    q <- sqrt(i*i - 2*i*f  + f*f)
    w <- exp(-(l + q) / 2)
    v <- exp(-(l - q) / 2)
    P <- array(dim=c(3, 3))
    P[1,1] <- (2*(v-w)*f + v*(q-l) + w*(q+l)) / (2*q)
    P[2,1] <- 0 
    P[3,1] <- 0 
    
    P[1,2] <- i*(v - w)/q
    P[2,2] <- -((2*f - l)*(v-w) - q*(v+w)) / (2*q)
    P[3,2] <- 0
    
    P[1,3] <- (-l*(v-w) - q*(v+w))/(2*q) + 1
    P[2,3] <- ((v-w)*(2*f - l) - q*(v+w))/(2*q) + 1
    P[3,3] <- 1
    P
}

state_probs_rem <- function(i, f, r)  {
    l <- i + r + f
    q <- sqrt(i*i + 2*i*r -  2*i*f  + r*r + 2*f*r + f*f)
    w <- exp(-(l + q) / 2)
    v <- exp(-(l - q) / 2)
    P <- array(dim=c(3, 3))

    P[1,1] <- (2*(v-w)*(f+r) + v*(q-l) + w*(q+l)) / (2*q)
    P[2,1] <- (v-w)*r/q
    P[3,1] <- 0 
    
    P[1,2] <- i*(v - w)/q
    P[2,2] <- -((2*(f+r) - l)*(v-w) - q*(v+w)) / (2*q)
    P[3,2] <- 0
    
    P[1,3] <- (-l*(v-w) - q*(v+w))/(2*q) + 1
    P[2,3] <- ((v-w)*(2*f - l) - q*(v+w))/(2*q) + 1
    P[3,3] <- 1
    P
}
