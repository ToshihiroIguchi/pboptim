\name{pboptim}
\alias{pboptim}
\title{Population Based Search}
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
method = c("DEO", "PSO", "GA"),
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
\item{method}{
A vector that specifies the optimization technique.
Default is c("DEO", "PSO", "GA").
SOMA can also be selected, but SOMA can not set the initial value.
}
\item{population}{The population size.}
\item{generation}{The maximum number of iterations.}
\item{trace}{Show progress of optimization calculation.}
\item{maximize}{For maximization, specify TRUE. Default is FALSE.}

\item{strategy}{
defines the Differential Evolution strategy used in the optimization procedure:
1: DE / rand / 1 / bin (classical strategy)
2: DE / local-to-best / 1 / bin (default)
3: DE / best / 1 / bin with jitter
4: DE / rand / 1 / bin with per-vector-dither
5: DE / rand / 1 / bin with per-generation-dither
6: DE / current-to-p-best / 1
any value not above: variation to DE / rand / 1 / bin: either-or-algorithm. Default strategy is currently 2.
}
\item{bs}{
if FALSE then every mutant will be tested against a member in the previous generation, and the best value will proceed into the next generation (this is standard trial vs. target selection). If TRUE then the old generation and NP mutants will be sorted by their associated objective function values, and the best NP vectors will proceed into the next generation (best of parent and child selection). Default is FALSE.
}
\item{CR}{
crossover probability from interval [0,1]. Default to 0.5.
}
\item{F}{
differential weighting factor from interval [0,2]. Default to 0.8.
}
\item{p}{
when strategy = 6, the top (100 * p)\% best solutions are used in the mutation. p must be defined in (0,1].
}
\item{c}{
c controls the speed of the crossover adaptation. Higher values of c give more weight to the current successful mutations. c must be defined in (0,1].
}

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
