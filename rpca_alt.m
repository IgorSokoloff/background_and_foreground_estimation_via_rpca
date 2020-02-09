
function [L, S, Rel_error, Time] = rpca_alt(A, r, alpha, params)

% A : A data matrix to be decomposed into a low-rank matrix L and a sparse
% matrix M. Unobserved entries are represented as zeros.
% r : Target rank
% alpha : An upper bound of max sparsity over the columns/rows of S
% params : parameters for the algorithm
%   .max_iter : Maximum number of iterations (default 30)
%   .tol : Desired Frobenius norm error (default 2e-4)
%%%%%%%%%%%
% Output:: L : Low rank matrix of rank r
%          S : Sparse matrix of sparsity alpha*m*n
%  Rel_error : ||A-L-S||_F/||A||_F

tic
% Default parameter settings

%%%% Read paramter settings
if isfield(params,'gamma')
    gamma = params.gamma;
end
if isfield(params,'max_iter')
    maxIter = params.max_iter;
end
if isfield(params,'thresh')
    thresh= params.thresh;
end

% Add Library paths
addpath PROPACK;
addpath bksvd-master;

%%%%%% Setting up %%%%%%%

Norm_of_A = norm(A,'fro');
[d1, d2] = size(A);

is_sparse  = issparse(A);
if is_sparse
    [I, J, A_vec] = find(A);
    n = length(A_vec);
    obs_ind = sub2ind([d1,d2], I, J);
    col = [0; find(diff(J)); n];
    p = n/d1/d2;
    if p>0.9
        is_sparse = 0;
        A = full(A);
    end
else
    p = 1;
end

%%%%% Initialization %%%%%%%
L_old = zeros(size(A));
S_old = zeros(size(A));
Rel_error = zeros(1, maxIter); converged = 0; %%(When showing converegcne of the algorithm)
Time = zeros(1, maxIter);
iter = 0;

%%%% Iteration %%%%%%%%%
%for i = 1:maxIter
while ~converged
    iter  = iter+1;

    %%% Projection on the convex set::: A=L+S
    L_temp = L_old/2-S_old/2+A/2;
    S_temp = S_old/2-L_old/2+A/2;

    %%%% Low rank Projection%%%%%%

    if p==1
        [U,Sigma,V] = bksvd(L_temp,r,1);
    else
        [U, Sigma, V] = lansvd(L_temp,r,'L');
    end

    L = U*Sigma*V';

    %%%% Sparse Projection %%%%%%

    S = Tproj(S_temp, gamma*p*alpha, gamma*p*alpha);

    Residual = A-L-S;
    Rel_error(1, iter) = norm(Residual,'fro')/Norm_of_A;
    Time (1, iter) = toc;
    L_old = L;
    S_old = S;

    %%%%  Convergence criteria %%%%
    if (Rel_error(iter) <= thresh ||iter>maxIter)
         converged = 1;
    end
    
end
toc