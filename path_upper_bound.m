% compute the bound for path
function ub = path_upper_bound(k, A, n, eindex)
    %% path specific upper bound
    % k is the number of edges in the path being added
    % A is the adjacent matrix
    sigs = 2*cos(pi*(1:k+1)/(k+2));
    sigs = sigs(sigs > 0);
    hk = length(sigs);
    esum = n*exp(eindex);
    % need top k eigenvalues of A
    lambda1 = abs(eigs(A,1));
    % eigs returns the largest *magnitude* eigenvalues. This transformastion
    % ensures we get the *algebraically* largest.
    lambs = eigs(A+lambda1*speye(n),hk)' - lambda1;
    esumnew = esum - sum(exp(lambs));
    esumnew = esumnew + sum(exp(lambs+sigs));
    % this is the final upper bound
    ub = log(esumnew/n);
end