//
//  bioinspired.cpp
//  proyecto00
//
//  Created by Jaime Francisco Aguayo on 13/07/18.
//  Copyright © 2018 Jaime Francisco Aguayo. All rights reserved.
//

#include "bioinspired.hpp"

/// Chicken Swarm Optimization
double cso(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration,
           int G, double FL, const char name[])
{
    gsl_rng *rng_ptr; // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    const double eps = 0.001;
    
    int rn = ceil(0.15 * n);    // Number of roosters
    int hn = ceil(0.7 * n);     // Number of heins
    int cn = n - rn - hn;       // Number of chicks
    int mn = ceil(0.2 * n);     // Numer of mother heins
    
    // Control of hierarchy
    int roosters[rn];
    int hens[hn];
    int chicks[cn];
    int mothers[mn];
    int chicksMum[cn];  // Save the mums of the chicks
    int hensLider[hn];// Save the hens' rooster lider
    
    // Setting chickens
    double ** agents = createMatrix(n, dimension);
    double ** pbest = createMatrix(n, dimension);
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    copyMatrix(agents, pbest, n, dimension);
    
    // Setting fitness
    double * fitness = createVector(n);
    for (int i=0; i<n; i++) {
        fitness[i] = function(agents[i]);
    }
    double * pfit = createVector(n);
    copyVector(fitness, pfit, n);
    
    // To save the best
    double * Pbest = createVector(dimension);
    double* i1; i1 = std::min_element(fitness, fitness+n);
    copyVector(agents[(i1-fitness)], Pbest, dimension);
    double * Gbest = createVector(dimension); copyVector(Pbest, Gbest, dimension);
    
    // Iterations
    for (int k=0; k<iteration; k++) {
        // If is time to update hierarchies
        if (k % G == 0) {
            update_relationship(n, function, agents, fitness, rn, hn, cn, mn, roosters, hens, chicks, mothers, chicksMum, hensLider);
        }
        
        // Search for food
        
        // Roosters
        for (int i=0; i<rn; i++) {
            int ri = roosters[i];
            // Pick a different rooster
            int r2 = roosters[(rand() % rn)];
            while (ri == r2) {
                r2 = roosters[(rand() % rn)];
            }
            
            // Compute variance
            double sigma;
            if (pfit[ri] <= pfit[r2])
                sigma = 1.0;
            else
                sigma = sqrt( exp((pfit[r2] - pfit[ri]) / (abs(pfit[ri]) + eps)) );
            
            // Updating position of the i-th rooster
            for (int d = 0; d<dimension; ++d)
                agents[ri][d] = pbest[ri][d]*(1.0 + gsl_ran_gaussian(rng_ptr, sigma));
            
        }
        
        // Hens
        for (int i=0; i<hn; i++) {
            int heni = hens[i];
            
            int henLider = hensLider[i];
            // Choose chicken to steal food
            int aux = (rand() % hn);
            int a = hens[aux]; int ar = hensLider[aux];
            int b = roosters[(rand() % rn)];
            int chicken = (rand()%2)? a : b ;
            ar = ((chicken == a)? ar:b);
            while (henLider == ar) {    // Not in the same group
                aux = (rand() % hn);
                a = hens[aux]; ar = hensLider[aux];
                b = roosters[(rand() % rn)];
                chicken = (rand()%2)? a : b ;
                ar = ((chicken == a)? ar:b);
            }
            
            // Calculate s1 and s2
            double s1 = exp((pfit[heni] - pfit[chicken]) / (abs(pfit[heni]) + eps));
            double s2 = exp(pfit[chicken] - pfit[heni]);
            if(isinf(s2)) s2 = DBL_MAX;
            
            double rand1;
            double rand2;
            // Update position of the ith hen
            for (int d = 0; d<dimension; ++d) {
                rand1 = gsl_ran_flat(rng_ptr, 0.0, 1.0);
                rand2 = gsl_ran_flat(rng_ptr, 0.0, 1.0);
                agents[heni][d] = pbest[heni][d]
                + s1 * rand1 * (pbest[henLider][d] - pbest[heni][d])
                + s2 * rand2 * (pbest[chicken][d] - pbest[heni][d]);
            }
        }
        // Chicks
        for (int i=0; i<cn; i++) {
            for (int d = 0; d<dimension; ++d)
                agents[chicks[i]][d] = pbest[chicks[i]][d] * FL * ( pbest[chicksMum[i]][d] - pbest[chicks[i]][d]);
        }
        
        // Checking constrains
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                if (agents[i][j] < lb[j])
                    agents[i][j] = lb[j];
                else if (agents[i][j] > ub[j])
                    agents[i][j] = ub[j];
            }
        }
        
        // Checking if there were progress
        for (int i=0; i<n; i++) {
            fitness[i] = function(agents[i]);
            if (fitness[i] < pfit[i]) {
                pfit[i] = fitness[i];
                copyVector(agents[i], pbest[i], dimension);
            }
            else
                fitness[i] = pfit[i];   // Si no fue mejor, ni le muevas
        }
        
        // Saving the best of the fitness
        i1 = std::min_element(fitness, fitness+n);
        copyVector(agents[(i1-fitness)], Pbest, dimension);
        if (function(Pbest) < function(Gbest))
            copyVector(Pbest, Gbest, dimension);
        if (name != NULL) {
            file<<k<<" "<<(*i1)<<"/n";
        }
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    freeMatrix(agents);
    freeMatrix(pbest);
    free(fitness);
    free(pfit);
    free(Pbest);
    free(Gbest);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}


void update_relationship(int n, double (*function)(double *), double ** agents, double * fitness, int rn, int hn, int cn, int mn, int roosters[], int hines[], int chicks[], int mothers[], int chicksMum[], int hinesLider[])
{
    // Numerate the chickens and sort them by fitness values
    std::vector< std::pair <double,int> > vect;
    for (int i=0; i<n; i++)
        vect.push_back( std::make_pair(fitness[i], i) );
    sort(vect.begin(), vect.end());
    
    for (int i=0; i<rn; i++)        // Saving roosters
        roosters[i] = vect[i].second;
    
    for (int i=rn; i<(rn+hn); i++)  // Saving hines
        hines[i-rn] = vect[i].second;
    
    for (int i=rn+hn; i<n; i++)     // Saving chicks
        chicks[i-rn-hn] = vect[i].second;
    
    // Shuffling hines to assign mothers
    std::random_shuffle(hines, hines+hn);
    for (int i=0; i<mn; i++)
        mothers[i] = hines[i];
    
    // Assign every chick to its mummy
    for (int i=0; i<cn; i++) {
        int rnd = rand() % mn;
        chicksMum[i] = mothers[rnd];
    }
    
    // Assign every chicken to a rooster to be its lider
    for (int i=0; i<hn; i++) {
        int rnd = rand() % rn;
        hinesLider[i] = roosters[rnd];
    }
    
}

// ------------------------ Cuckoo Search Algorithm ------------------------------- //

double csa(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, int nest, double pa, const char name[])
{
    gsl_rng *rng_ptr; // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    double ** agents = createMatrix(n, dimension);
    double ** nests = createMatrix(nest, dimension);
    double * step = createVector(dimension);
    double * stepsize = createVector(dimension);
    
    // To create the Lévy flights
    double beta = 3.0 / 2.0;
    double sigma = (tgamma(1.0 + beta) *
                    sin(3.141592 * beta / 2.0) / (tgamma((1.0 + beta) / 2.0) * beta * pow(2.0, ((beta - 1.0) / 2.0)) ) );
    sigma = pow(sigma, (1.0 / beta));
    
    for (int i=0; i<dimension; i++) {
        double u = gsl_ran_gaussian(rng_ptr, 1.0) * sigma;
        double v = gsl_ran_gaussian(rng_ptr, 1.0);
        step[i] = u / pow (abs(v), (1.0 / beta));
    }
    
    // Setting agents
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    
    // Setting nests
    for (int i=0; i<nest; i++) {
        for (int j=0; j<dimension; j++) {
            nests[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
        }
    }
    
    // Saving fitness of nests
    double * fitness = createVector(nest);
    for (int i=0; i<nest; i++) {
        fitness[i] = function(nests[i]);
    }
    
    // To save the best
    double * Pbest = createVector(dimension);
    double* i1; i1 = std::min_element(fitness, fitness+nest);
    copyVector(nests[(i1-fitness)], Pbest, dimension);
    double * Gbest = createVector(dimension); copyVector(Pbest, Gbest, dimension);  // The global best
    
    for (int t = 0; t<iteration; t++) {                         // Iterate
        for (int i=0; i<n; i++) {                               // For every agent
            int val = rand() % nest;
            double aux;// Take a nest
            if ((aux=function(agents[i])) < fitness[val]) {     // If agent fitness better than nest's
                copyVector(agents[i], nests[val], dimension);   // agent -> nest
                fitness[val] = aux;                             // Update fitness
            }
        }
        
        // Delete nests with worst fitness with prob. pa
        delete_worst(function, nests, agents, fitness, n, nest, dimension, pa, lb, ub, rng_ptr);
        
        // Checking constrains for nests
        for (int i=0; i<nest; i++) {
            for (int j=0; j<dimension; j++) {
                if (nests[i][j] < lb[j])
                    nests[i][j] = lb[j];
                else if (nests[i][j] > ub[j])
                    nests[i][j] = ub[j];
            }
        }
        
        // Update agents (Lévy flights)
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                stepsize[j] = 0.2 * step[j] * (agents[i][j] - Pbest[j]);
                agents[i][j] += stepsize[j] * gsl_ran_gaussian(rng_ptr, 1.0);
            }
        }
        
        // Checking constrains for agents
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                if (agents[i][j] < lb[j])
                    agents[i][j] = lb[j];
                else if (agents[i][j] > ub[j])
                    agents[i][j] = ub[j];
            }
        }
        
        // To save the bests nests
        for (int i=0; i<nest; i++) {
            fitness[i] = function(nests[i]);
        }
        Pbest = createVector(dimension);
        i1 = std::min_element(fitness, fitness+nest);
        copyVector(nests[(i1-fitness)], Pbest, dimension);
        
        if (function(Pbest) < function(Gbest))
            copyVector(Pbest, Gbest, dimension);
        if (name != NULL) {
            file<<t<<" "<<(*i1)<<"/n";
        }
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    freeMatrix(agents);
    freeMatrix(nests);
    free(step);
    free(stepsize);
    free(fitness);
    free(Pbest);
    free(Gbest);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}

void delete_worst(double (*function)(double *), double ** nests, double ** agents, double *fitness, int n, int nest, int dim, double pa, double * lb, double * ub, gsl_rng *rng_ptr)
{
    // Numerate the nest and sort them by fitness values
    std::vector< std::pair <double,int> > vect1;
    for (int i=0; i<nest; i++)
        vect1.push_back( std::make_pair(fitness[i], i) );
    sort(vect1.begin(), vect1.end());
    
    // Numerate the cuckoos and sort them by fitness values
    std::vector< std::pair <double,int> > vect2;
    for (int i=0; i<n; i++)
        vect2.push_back( std::make_pair(function(agents[i]), i) );
    sort(vect2.begin(), vect2.end());
    
    int nworst = nest / 2;
    // Delete worst fitness with probability pa
    for (int i=nest-1; i>=nworst; i--) {
        if (gsl_ran_flat(rng_ptr, 0.0, 1.0) < pa)
            for (int j=0; j<dim; j++)
                nests[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    
    int mworst;
    if (nest > n)
        mworst = n;
    else
        mworst = nest;
    
    // Keeping the best nests
    for (int i=0; i<mworst; i++) {
        if (vect1[i].first < vect2[i].first)
            copyVector(nests[vect1[i].second], agents[vect2[i].second], dim);
    }
    
}

/// Algoritmo revisado
double abc2(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, const char name[])
{
    gsl_rng *rng_ptr; // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    n = (n + n % 2);
    
    double ** agents = createMatrix(n, dimension);
    double * fitness = createVector(n);
    int trials[n];
    double * probabilities = createVector(n);
    double * Gbest = createVector(dimension);
    
    
    // Bees at random position
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    
    // Setting fitness
    for (int i=0; i<n; i++) {
        fitness[i] = function(agents[i]);
    }
    
    // Saving the bee with best fitness
    double* i1; i1 = std::min_element(fitness, fitness+n);
    copyVector(agents[(i1-fitness)], Gbest, dimension);
    
    double beta = 0.0;
    for (int t=0; t<iteration; t++) {
        // Employees phase
        for (int i=0; i<n; i++)
            send_employees(n, dimension, i, agents, lb, ub, function, fitness, trials);
        
        // Calculate probabilities
        i1 = std::max_element(fitness, fitness+n);
        for (int i=0; i < n; i++)
            probabilities[i] = 1.0 - (0.9 * fitness[i] / (*i1) + 0.1);
        
        // Onlookers phase
        double phi = gsl_ran_flat(rng_ptr, 0.0, 1.0);
        i1 = std::max_element(probabilities, probabilities+n);
        beta += phi * (*i1);
        beta = abs(beta - (*i1)*(int)(beta/(*i1)));
        
        // Sending onlooker (s)
        double sum = 0.0;
        int index = 0;
        for (int i=0; i < n; i++) {
            sum += probabilities[i];
            if (beta < sum)
                index = i;
        }
        send_employees(n, dimension, index, agents, lb, ub, function, fitness, trials);
        
        // Scout phase
        for (int i=0; i < n; i++) {
            int maxTrials = (int)(0.6 * n * dimension);
            if (trials[i] > maxTrials) {
                for (int j=0; j<dimension; j++)
                    agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
                fitness[i] = function(agents[i]);
                trials[i] = 0;
            }
        }
        
        i1 = std::min_element(fitness, fitness+n);
        if (fitness[(i1-fitness)] < function(Gbest))
            copyVector(agents[(i1-fitness)], Gbest, dimension);
        if (name != NULL) {
            file<<t<<" "<<(*i1)<<"/n";
        }
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    
    freeMatrix(agents);
    free(probabilities);
    free(fitness);
    free(Gbest);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}


// Crea nuevas abejas que irán tras de las mejores
void send_employees(int n, int dim, int i, double ** agents, double * lb, double * ub, double (*function)(double *), double * fitness, int trials[])
{
    gsl_rng *rng_ptr;                           // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    
    double mutant[dim];
    
    int d = rand() % dim;
    int j = i;
    while (j == i)
        j = rand() % n;
    // Mutate
    mutant[d] = agents[i][d] + (gsl_ran_flat(rng_ptr, lb[d], ub[d]) - 0.5) * 2.0 * (agents[i][d] - agents[j][d]);
    if(mutant[d] < lb[d])
        mutant[d] = lb[d];
    else if (mutant[d] > ub[d])
        mutant[d] = ub[d];
    
    if (function(mutant) < function(agents[i])) {
        copyVector(mutant, agents[i], dim);
        fitness[i] = function(agents[i]);
        trials[i] = 0;
    }
    else
        trials[i]++;
    
    gsl_rng_free(rng_ptr);
}

// Gray Wolfe Optimization
double gwo(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, const char name[])
{
    gsl_rng *rng_ptr; // pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    double ** agents = createMatrix(n, dimension);
    double * Gbest = createVector(dimension);
    
    int alpha, beta, delta;
    
    // Wolfes at random position
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    
    // Ranking the wolfes by fitness
    std::vector< std::pair <double,int> > ranking;
    for (int i=0; i<n; i++)
        ranking.push_back( std::make_pair(function(agents[i]), i) );
    sort(ranking.begin(), ranking.end());
    
    // Setting alpha, beta, delta
    alpha = ranking[0].second;
    beta = ranking[1].second;
    delta = ranking[2].second;
    
    // Setting the best
    copyVector(agents[alpha], Gbest, dimension);
    
    for (int t=0; t < iteration; t++) {
        double a = 2.0 - (2.0 * t) / iteration;     // a decreases linearly fron 2 to 0
        double r1, r2, A1, A2, A3, C1, C2, C3, D_alpha, D_beta, D_delta, X1, X2, X3;
        
        // Update positions
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                r1 = gsl_ran_flat(rng_ptr, 0, 1);   // r1 is a random number in [0,1]
                r2 = gsl_ran_flat(rng_ptr, 0, 1);   // r2 is a random number in [0,1]
                A1 = 2.0 * a * r1 - a;              // Equation (3.3)
                C1 = 2.0 * r2;                      // Equation (3.4)
                D_alpha = abs(C1*agents[alpha][j] - agents[i][j]);  // Equation (3.5)
                X1 = agents[alpha][j] - A1*D_alpha; // Equation (3.6)
                
                r1 = gsl_ran_flat(rng_ptr, 0, 1);   // r1 is a random number in [0,1]
                r2 = gsl_ran_flat(rng_ptr, 0, 1);   // r2 is a random number in [0,1]
                A2 = 2.0 * a * r1 - a;              // Equation (3.3)
                C2 = 2.0 * r2;                      // Equation (3.4)
                D_beta = abs(C2*agents[beta][j] - agents[i][j]);    // Equation (3.5)
                X2 = agents[beta][j] - A2*D_beta;   // Equation (3.6)
                
                r1 = gsl_ran_flat(rng_ptr, 0, 1);   // r1 is a random number in [0,1]
                r2 = gsl_ran_flat(rng_ptr, 0, 1);   // r2 is a random number in [0,1]
                A3 = 2.0 * a * r1 - a;              // Equation (3.3)
                C3 = 2.0 * r2;                      // Equation (3.4)
                D_delta = abs(C3*agents[delta][j] - agents[i][j]);  // Equation (3.5)
                X3 = agents[delta][j] - A3*D_delta; // Equation (3.6)
                
                agents[i][j] = (X1+X2+X3) / 3.0;    // Equation (3.7)
            }
        }
        // Checking constrains for agents
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                if (agents[i][j] < lb[j])
                    agents[i][j] = lb[j];
                else if (agents[i][j] > ub[j])
                    agents[i][j] = ub[j];
            }
        }
        // Ranking the new positions
        ranking.clear();
        ranking.shrink_to_fit();
        for (int i=0; i<n; i++)
            ranking.push_back( std::make_pair(function(agents[i]), i) );
        sort(ranking.begin(), ranking.end());
        
        // Setting alpha, beta, delta
        alpha = ranking[0].second;
        beta = ranking[1].second;
        delta = ranking[2].second;
        
        double best = function(Gbest);
        if (ranking[0].first < best)
            copyVector(agents[alpha], Gbest, dimension);
        if (name != NULL) {
            file<<t<<" "<<best<<"/n";
        }
    }
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    freeMatrix(agents);
    free(Gbest);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}

/// Whale Swarm Optimization
double wso(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, double ro0, double eta, const char name[])
{
    gsl_rng *rng_ptr; // Pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    double ** agents = createMatrix(n, dimension);
    double ** new_agents = createMatrix(n, dimension);
    double * Gbest = createVector(dimension);
    
    // Whales at random position
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
    }
    
    // Ranking the whales by fitness
    std::vector< std::pair <double,int> > ranking;
    for (int i=0; i<n; i++)
        ranking.push_back( std::make_pair(function(agents[i]), i) );
    sort(ranking.begin(), ranking.end());
    copyVector(agents[ranking[0].second], Gbest, dimension);        // Saving the best
    copyMatrix(agents, new_agents, n, dimension);
    
    // Iterations
    for (int t=0; t<iteration; t++) {
        double dist;
        // For every Whale
        for (int i=0; i<n; i++) {
            int bnw = better_nearest_whale(i, ranking, agents, n, dimension, dist); // Search for the better n nearest whale
            if (bnw != -1) {                                        // If such whale exist
                for (int j=0; j<dimension; j++)
                    new_agents[i][j] = agents[i][j] + gsl_ran_flat(rng_ptr, 0.0, ro0 * exp(-eta*dist/20.0))
                    * (agents[bnw][j] - agents[i][j]); // Modify its pos.
            }
        }
        
        // Checking constrains for new whales
        for (int i=0; i<n; i++) {
            for (int j=0; j<dimension; j++) {
                if (new_agents[i][j] < lb[j])
                    new_agents[i][j] = lb[j];
                else if (new_agents[i][j] > ub[j])
                    new_agents[i][j] = ub[j];
            }
        }
        copyMatrix(new_agents, agents, n, dimension);               // Copy the new whales
        
        // Ranking the new positions
        ranking.clear();
        ranking.shrink_to_fit();
        for (int i=0; i<n; i++)
            ranking.push_back( std::make_pair(function(agents[i]), i) );
        sort(ranking.begin(), ranking.end());
        
        // Saving the best whale
        double best = function(Gbest);
        if (ranking[0].first < best)
            copyVector(agents[ranking[0].second], Gbest, dimension);
        
        if (name != NULL) {
            file<<t<<" "<<best<<"/n";
        }
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    gsl_rng_free(rng_ptr);
    freeMatrix(new_agents);
    freeMatrix(agents);
    free(Gbest);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}


int better_nearest_whale(int i, std::vector< std::pair <double,int> > ranking, double ** agents, int n, int dim, double &dist)
{
    int bnw = -1;
    double aux = DBL_MAX;
    dist = DBL_MAX;
    
    // For every whale better than ranking[i].second
    for (int j=0; j<i; j++) {
        aux = vectorDifNorm(agents[ranking[j].second], agents[ranking[i].second], dim);
        // If it is the closest so far, keep it
        if (aux < dist) {
            bnw = ranking[j].second;
            dist = aux;
        }
    }
    return bnw;
}


/// Monarch Butterfly Optimization
double mbo(int n, double (*function)(double *), double *x, double * lb, double * ub, int dimension, int iteration, double bar, double p, double peri, const char name[])
{
    gsl_rng *rng_ptr; // Pointer to random number generator (rng)
    rng_ptr = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng_ptr, 121);
    std::srand ( unsigned(std::time(0)) );
    std::ofstream file;
    if (name != NULL) {
        file = std::ofstream(name);
    }
    
    int nLand1 = ceil(p*n);             // NP1 in paper
    int nLand2 = n - nLand1;            // NP2 in paper
    int land1[nLand1];                  // To save the index of mb in land 1
    int land2[nLand2];                  // To save the index of mb in land 2
    
    double ** agents = createMatrix(n, dimension);
    double ** tempAgents = createMatrix(n, dimension);
    double * Gbest = createVector(dimension);
    double *fitness = createVector(n);
    
    // Butterflies at random position
    for (int i=0; i<n; i++) {
        for (int j=0; j<dimension; j++)
            agents[i][j] = gsl_ran_flat(rng_ptr, lb[j], ub[j]);
        fitness[i] = function(agents[i]);                       // Setting fitness
    }
    
    // Ranking the butterfly by fitness
    std::vector< std::pair <double,int> > ranking;
    for (int i=0; i<n; i++)
        ranking.push_back( std::make_pair(fitness[i], i) );
    sort(ranking.begin(), ranking.end());
    copyVector(agents[ranking[0].second], Gbest, dimension);    // Saving the best
    
    for (int t=0; t<iteration; t++) {
        copyMatrix(agents, tempAgents, n, dimension);
        /* Divide the whole population into Population1(Land1) and Population2(Land2)
         according to their fitness (Better butterflies in land 1). */
        for (int i = 0; i < n; ++i) {
            if (i < nLand1)
                land1[i] = ranking[i].second;
            else
                land2[i-nLand1] = ranking[i].second;
        }
        
        // Migration operator (algorithm 1)
        for (int i = 0; i < nLand1; i++) {
            for (int d = 0; d < dimension; d++) {
                if (gsl_ran_flat(rng_ptr, 0.0, 1.0) * peri <= p) {
                    int r1 = rand() % nLand1;                       // Pick a butterfly in land 1
                    agents[land1[i]][d] = tempAgents[land1[r1]][d]; // Eq 1
                }
                else {
                    int r2 = rand() % nLand2;                       // Pick a butterfly in land 2
                    agents[land1[i]][d] = tempAgents[land2[r2]][d]; // Eq 3
                }
            }
        }
        
        // Butterfly adjusting operator (algorithm 2)
        for (int i = 0; i < nLand2; i++)
        {
            double scale = 1.0 / pow(t+1, 2);                       // Smaller step for local walk
            
            for (int d = 0; d<dimension; d++) {
                double alpha = gsl_ran_flat(rng_ptr, 0.0, 2.0), beta = gsl_ran_flat(rng_ptr, -1.0, 1.0);
                double dx = gsl_ran_levy_skew(rng_ptr, 1.0, alpha, beta);   // Levy flight
                
                if (gsl_ran_flat(rng_ptr, 0.0, 1.0) >= p)
                    agents[land2[i]][d] = tempAgents[ranking[0].second][d];
                else {
                    int r3 = rand() % nLand2;                       // Pick a butterfly in land 2
                    agents[land2[i]][d] = tempAgents[land2[r3]][d]; // Eq. 4
                    if (gsl_ran_flat(rng_ptr, 0.0, 1.0) > bar)
                        agents[land2[i]][d] += scale*(dx - 0.5);// Eq. 6
                }
            }
            
        }
        
        // Checking constrains, evaluate fitness and ranking the new butterflies
        ranking.clear();
        ranking.shrink_to_fit();
        for (int i=0; i<n; i++) {                               // Checking constrains
            for (int j=0; j<dimension; j++) {
                if (agents[i][j] < lb[j])
                    agents[i][j] = lb[j];
                else if (agents[i][j] > ub[j])
                    agents[i][j] = ub[j];
            }
            fitness[i] = function(agents[i]);                   // Eval. fitness
            ranking.push_back( std::make_pair(fitness[i], i) ); // To ranking them
        }
        sort(ranking.begin(), ranking.end());                   // Sort by fitness
        
        // Save the best
        double best = function(Gbest);
        if (ranking[0].first < best)
            copyVector(agents[ranking[0].second], Gbest, dimension);
        if (name != NULL) {
            file<<t<<" "<<best<<"/n";
        }
    }
    
    // Answer
    copyVector(Gbest, x, dimension);
    double bestValue = function(Gbest);
    
    freeMatrix(agents);
    freeMatrix(tempAgents);
    free(fitness);
    free(Gbest);
    gsl_rng_free(rng_ptr);
    
    if (name != NULL) {
        file.close();
    }
    return bestValue;
}
