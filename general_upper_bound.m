function ub = general_upper_bound(k, A, n, natconn)
    %% gereral upper bound for arbitrary edge additions
    esum = exp(natconn);
    % need top eigenvalue of A
    lambda1 = abs(eigs(A,1));
    % need top 2k eigenvalues of A
    lambda1 = abs(eigs(A,1));
    lambs = eigs(A+lambda1*speye(n),2*k)' - lambda1;
    lambssum = sum(exp(lambs));
    factor = exp(sqrt(2*k)) + (2*k-1);
    % this is the final upper bound
    ub = log(esum - lambssum/n + factor*exp(lambda1)/n);
end