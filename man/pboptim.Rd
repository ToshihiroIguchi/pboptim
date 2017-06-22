\name{pboptim}
\alias{pboptim}
\title{Population Based Search.}
\description{
Pboptim is a wrapper for population based search.
Algorithms "Differential Evolution Optimization",
"Particle Swarm Optimization",
"Genetic Algorithms" and
"Self-Organising Migrating Algorithm" can be used.
"pboptim" depends on "DEoptim", "pso", "GA", and "soma" packages.
}
\usage{
pboptim(fn ,lower, upper, initialpar = NULL,
method = c("DEO", "PSO", "GA", "SOMA"),
population = 20, generation = 10,
trace = TRUE, maximize = FALSE,
DEOcontrol = list(strategy = 2, bs = FALSE, CR = 0.5, F = 0.8, p = 0.2, c = 0),
PSOcontrol = list(s = floor(10+2*sqrt(length(pso_par))),
k = 3, w = 0.7213, c.p = 1.193, type = "SPSO2007"),
GAcontrol = list(pcrossover = 0.8,pmutation = 0.1,
elitism = base::max(1, round(population * 0.05))),
SOMAcontrol = list(pathLength = 3, stepLength = 0.11, perturbationChance = 0.1))

}

\arguments{
\item{fn}{A function to be minimized (or maximized).}
\item{lower}{Lower bounds on the variables.}
\item{upper}{Upper bounds on the variables.}
\item{initialpar}{A vector representing the initial value. Default is NULL.}
\item{method}{A vector that specifies the optimization technique.}
\item{population}{The population size.}
\item{generation}{The maximum number of iterations.}
\item{trace}{Show progress of optimization calculation.}
\item{maximize}{For maximization, specify TRUE. Default is FALSE.}

}

\examples{

#Population Based Search
result <- pboptim(function(x){x[1]^2+x[2]^2},
                method = c("DEO", "PSO", "GA", "SOMA"),
                initialpar = c(1, 1),
                lower = c(-2,-2), upper = c(2,2),
                population = 20, generation = 100,
                maximize = FALSE)

#View generation vs value
plot(result)

#summary
summary(result)
}

\references{
Cortez, Paulo. Modern optimization with R. Springer, 2014.
}
