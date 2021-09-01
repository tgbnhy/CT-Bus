function [trace_est, var]= deltashift_hutchpp(A1, A2, prev_tr, prev_var, num_queries, iter)
    %Use parameter-free deltashift with hutch++ to estimate trace of exp(A2)
    % A1: Previous matrix, A2: Current matrix in the sequence
    % prev_tr: Previous trace estimate (of A1)
    % prev_var: Previous variance estimate
    % num_queries: Number of mat-vec to be used
    % iter: Lanczos iterations
    
    %   Hutch++ code adapted from https://github.com/RaphaelArkadyMeyerNYU/HutchPlusPlus
    
    % Calculate which matrices get how many queries, and generate random sign matrices
    S = 2*randi(2,size(A1,1),ceil(num_queries/3))-3;
    G = 2*randi(2,size(A1,1),floor(num_queries/3))-3;
    l_s = size(S,2);
    l_g = size(G,2);
    
    %Trace and variance of Delta = A2-A1
    AS1 = lanczos(A1, S, @exp, iter);
    AS2 = lanczos(A2, S, @exp, iter);
    
    Delta_S = AS2 - AS1;
    [Q_Delta,~] = qr(Delta_S,0);
    G_Delta = G - Q_Delta*(Q_Delta'*G);
    
    AQ1 = lanczos(A1, Q_Delta, @exp, iter);
    AQ2 = lanczos(A2, Q_Delta, @exp, iter);
    
    Delta_Q = AQ2 - AQ1;
    AG1 = lanczos(A1, G_Delta, @exp, iter);
    AG2 = lanczos(A2, G_Delta, @exp, iter);
    
    Delta_G = AG2 - AG1;

    %Trace
    trace_delta = trace(Q_Delta'*Delta_Q) + (1/size(G_Delta,2))*trace(G_Delta'*Delta_G);
    %Variance
    K_Delta = Delta_S - (Delta_Q * Q_Delta' * S);
    N = trace(K_Delta' * K_Delta)/l_s;

    %Trace and variance of A2
    [Q2,~] = qr(AS2,0);
    G2 = G - Q2*(Q2'*G);
    AQ2_ = lanczos(A2, Q2, @exp, iter);
    AG2_ = lanczos(A2, G2, @exp, iter);
    
    trace_A2 = trace(Q2'*AQ2_) + (1/size(G2,2))*trace(G2'*AG2_);
    
    K_2 = AS2 - (AQ2_ * Q2' * S);
    M = trace(K_2' * K_2)/l_s;
    %Var_2 = (4/l) * M;

    %Optimal gamma
    gamma = (4*N + l_g*prev_var)/(4*M + 4*N + l_g*prev_var);
    
    %Estimates
    trace_est = (gamma * trace_A2) + (1-gamma)*(trace_delta + prev_tr);
    var = (4*gamma^2/l_g)*M + (4*(1-gamma)^2/l_g)*N + (1-gamma)^2 * prev_var;
    
    
end  