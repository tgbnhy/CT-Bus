function [trace_est, var]= deltashift_hutch(A1, A2, prev_tr, prev_var, num_queries, iter)
    %Use parameter-free deltashift to estimate trace of exp(A2)
    % A1: Previous matrix, A2: Current matrix in the sequence
    % prev_tr: Previous trace estimate (of A1)
    % prev_var: Previous variance estimate
    % num_queries: Number of mat-vec to be used
    % iter: Lanczos iterations
    
    
    n = size(A1,1);
    G = 2*randi(2,n,num_queries)-3;
    
    %Approximate exp(A) * G
    z = lanczos(A1,G,@exp,iter);
    w = lanczos(A2,G,@exp,iter);
    
    
    %Optimal gamma
    N = (1/num_queries) * trace(z' * z);
    M = (1/num_queries) * trace(w' * w);
    C = (1/num_queries) * trace(w' * z);
    gamma = 1 - ( 2*C / (num_queries*prev_var + 2*N) );
    
    trace_est = (1-gamma)* prev_tr + ( (1/num_queries) * trace(G' * (w - (1-gamma)*z) ) );
    var = (1-gamma)^2 * prev_var + (2/num_queries) * (N + (1-gamma)^2 * M - 2*(1-gamma)*C);
