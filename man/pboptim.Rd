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
trace = TRUE, maximize = FALSE)
}

\arguments{
\item{fn}{A function to be minimized (or maximized).}
\item{lower}{Lower bounds on the variables.}
\item{upper}{Upper bounds on the variables.}
\item{initialpar}{A vector representing the initial value. Default is NULL.}
\item{method}{A vector that specifies the optimization technique.}
\item{population}{The population size.}
\item{generation}{The maximum number of iterations.}
\item{trace}{}
\item{maximize}{For maximization, specify TRUE. Default is FALSE.}

}

\examples{
result <- pboptim(function(x){x[1]^2+x[2]^2},
                method = c("DEO", "PSO", "GA", "SOMA"),
                initialpar = c(1, 1),
                lower = c(-2,-2), upper = c(2,2),
                population = 20, generation = 100,
                maximize = FALSE)

plot(result)
summary(result)
}
