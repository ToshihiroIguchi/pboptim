## pboptim package
Population Based Search

### Description
Pboptim is an R package that performs multiple population based search.
Algorithms Differential Evolution Optimization,
Particle Swarm Optimization,
Genetic Algorithms and
Random Search can be used.
"pboptim" depends on "DEoptim", "pso", and "GA" packages.

### Installation
You can install from R console.

If devtools is not installed on your PC, install devtools with Internet connection.

    install.packages("devtools")

Install from GitHub using devtools.
    
    library(devtools)
    install_github("ToshihiroIguchi/pboptim")

Load the pboptim package and attach it.

    library(pboptim)


### Examples

First define the function you want to optimize.

    fn <- function(x){ x[1]^2 + x[2]^2}

Perform optimization.
    
    result <- pboptim(fn,
    method = c("DEO", "PSO", "GA"),
    initialpar = c(1, 1),
    lower = c(-2,-2), upper = c(2,2),
    population = 20, generation = 100,
    maximize = FALSE)


View generation vs value.

    plot(result)

summary

    summary(result)


### License 
MIT

### Auther
Toshihiro Iguchi

