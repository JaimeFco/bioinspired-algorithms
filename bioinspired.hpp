//
//  bioinspired.hpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 13/07/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#ifndef bioinspired_hpp
#define bioinspired_hpp

#include <stdio.h>
#include "metodos_numericos.hpp"
#include <vector>



// Chicken Swarm Optimization
//  n: number of agents
//  function: test function
//  lb: lower limits for plot axes
//  ub: upper limits for plot axes
//  dimension: space dimension
//  iteration: number of iterations
//  G: Number of rosters
//  FL: percent of chicks
double cso(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, int G=5, double FL=0.5, const char name[] = NULL);

void update_relationship(int n, double (*function)(double *), double ** agents, double * fitness, int rn, int hn, int cn, int mn, int roosters[], int hines[], int chicks[], int mothers[], int chicksMum[], int hinesLider[]);


///   Cuckoo Search Algorithm
//    n: number of agents
//    function: test function
//    lb: lower limits for plot axes
//    ub: upper limits for plot axes
//    dimension: space dimension
//    iteration: number of iterations
//    nest: number of nests (default value is 70)
//    pa: probability of cuckoo's egg detection (default value is 0.25)
double csa(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, int nest=50, double pa=0.25, const char name[] = NULL);

void delete_worst(double (*function)(double *), double ** nests, double ** agents, double *fitness, int n, int nest, int dim, double pa, double * lb, double * ub, gsl_rng *rng_ptr);



/// Artificial Bee Colony
//  n: number of agents
//  function: test function
//  lb: lower limits for plot axes
//  ub: upper limits for plot axes
//  dimension: space dimension
//  iteration: number of iterations
double abc2(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, const char name[] = NULL);

void send_employees(int n, int dim, int i, double ** agents, double * lb, double * ub, double (*function)(double *), double * fitness, int trials[]);

/// Gray Wolfe Optimization
//  n: number of agents
//  function: test function
//  lb: lower limits for plot axes
//  ub: upper limits for plot axes
//  dimension: space dimension
//  iteration: number of iterations
double gwo(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, const char name[] = NULL);

/// Whale Swarm Optimization
//  n: number of agents
//  function: test function
//  lb: lower limits for plot axes
//  ub: upper limits for plot axes
//  dimension: space dimension
//  iteration: the number of iterations
//  ro0: intensity of ultrasound at the origin of source (default value is 2)
//  eta: probability of message distortion at large distances (default value is 0.005)
double wso(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, double ro0 = 2.0, double eta = 0.05, const char name[] = NULL);

int better_nearest_whale(int i, std::vector< std::pair <double,int> > ranking, double ** agents, int n, int dim, double &dist);

/// Monarch Butterfly Optimization
//  n: number of agents
//  function: test function
//  lb: lower limits for plot axes
//  ub: upper limits for plot axes
//  dimension: space dimension
//  iteration: the number of iterations
//  bar: Butterfly Adjusting Rate (default value is 0.41667)
//  p: Partition value (default value is 0.41667)
//  peri: Period value (default value is 2)
double mbo(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, double bar = 0.4166, double p = 0.41667, double peri = 1.2, const char name[] = NULL);

#endif /* bioinspired_hpp */
