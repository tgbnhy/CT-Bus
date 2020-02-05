%compute the natural connectivity using stochastic lanczos method.
function conn = natural_connectivity(A, n, b, reps, iter)
    % Allocate space for Krylov subspace and tridiagonal matrix.
    K1 = b;
    K2 = zeros(n,reps);
    K3 = zeros(n,reps);
    beta = zeros(iter,reps);
    alpha = zeros(iter,reps);
    
    K2 = A*K1;%this can be optimized
    alpha(1,:) = dot(K2,K1);
    for i = 1:iter-1
        K2 = K2 - bsxfun(@times,K1,alpha(i,:));
        beta(i+1,:) = sqrt(dot(K2,K2));
        K2 = bsxfun(@times,K2,1./beta(i+1,:));%functions: 
        K3 = A*K2 - bsxfun(@times,K1,beta(i+1,:));
        alpha(i+1,:) = dot(K3,K2);
        K1 = K2; K2 = K3;
    end

    ests = zeros(1,reps);
    for z = 1:reps
        T = spdiags([[beta(2:iter,z)',0]' alpha(:,z) beta(:,z)], -1:1, iter, iter);%Return a sparse matrix from diagonals.
        [U,S] = eig(full(T)); %returns a vector of the eigenvalues of matrix A.
        ests(z) = U(1,:)*diag(exp(diag(S)))*U(1,:)';
    end            
    conn = log(mean(ests));
end