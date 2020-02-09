
%thresh - eps
function [L, S, Rel_error, Time] = rpca_alf(D, r, params)
tic
Norm_of_D = norm(D,'fro');

[d1, d2] = size(D);

%%%% Read paramter settings
if isfield(params,'lambda')
    lambda = params.lambda;
end
if isfield(params,'beta')
    beta = params.beta;
end
if isfield(params,'gamma')
    gamma = params.gamma;
end
if isfield(params,'max_iter')
    maxIter = params.max_iter;
end
if isfield(params,'thresh')
    thresh= params.thresh;
end
if isfield(params,'rho_0')
    rho_old= params.rho_0;
end
if isfield(params,'rho_max')
    rho_max= params.rho_max;
end
if isfield(params,'sur_kind')
    sur_kind= params.sur_kind;
end

% Add Library paths
addpath PROPACK;

%%%%% Initialization %%%%%%%
L_old = zeros(size(D));
Lambda_old = zeros(size(D));

V_old = eye([d2, r]);
S_old = zeros(size(D));

Rel_error = zeros(1, maxIter);
Time = zeros(1, maxIter);
converged = 0; %%(When showing converegcne of the algorithm)
iter = 0;


%%%% Iteration %%%%%%%%%
%for i = 1:maxIter
while ~converged
    iter  = iter+1;
    
    T = D + Lambda_old/rho_old;
    TS_dif = T - S_old;
    
    [A, ~, B] = lansvd(TS_dif*V_old, r, 'L');
    %[A, ~, B] = svd( TS_dif*V_old );
    
    U_temp = A*B';
    
    [A_v, Sigma_v, B_v] = lansvd(V_old, r, 'L');
    d_v = diag(Sigma_v);
    
    l = get_l(sur_kind, d_v, gamma); %equation (18) in the paper
    
    %disp ([size(A_v), size(l), size(B_v)]);
    
    assert (isvector(l));
    subgrad_V = A_v*diag(l)*B_v';
    
    V_temp = TS_dif'*U_temp - lambda*subgrad_V./rho_old ;
    L_temp = U_temp*V_temp'; 
    
    W = T - L_temp; 
    S_temp = max(abs(W)-1/rho_old, 0 ).*sign (W);
    
    Residual = D - L_temp - S_temp;
    
    Lambda_temp = Lambda_old + rho_old.* Residual;
    rho_temp = min(beta*rho_old, rho_max);
     
    Rel_error(1, iter) = norm(Residual, 'fro')/Norm_of_D;
    Time (1, iter) = toc;
    %disp ([iter, Rel_error(iter), rho_temp]);
    V_old = V_temp;
    S_old = S_temp;
    L_old = L_temp;
    Lambda_old = Lambda_temp;
    rho_old = rho_temp;
    
%%%%  Convergence criteria %%%%
    if (Rel_error(iter) <= thresh ||iter>maxIter)
         converged = 1;
    end
end
    L = L_temp;
    S = S_temp;
    
toc
end