function [X, W, mle_value,  fvals , dvals ] = dpglasso (S,maxiter,tol,ret_obj,Ydata,Xdata,mu,A,B,O,r,k_B)
% dpglasso Solve the graphical lasso problem using the DP-glasso algorithm .
%
% S is the  sample covariance matrix.
% r>0 is the regularization parameter .
% Returns precision matrix X and covariance matrix W.
% Returns primal and dual objective function values if retobj is true .
%
% File : dpglasso.m
% By default do not return objective function values .

if nargin < 5
    ret_obj = false;
end


% Get the dimension of the problem .
p = size(S, 1);


% Initialize the precision and covariance matrices.
% % dia = diag(S) + r;
% dia = diag(S);
% X = spdiags(1./dia, 0, p, p);
% W = diag(dia); % Covariance estimate will not be sparse in general.


% Dual feasible point.
W = S;
X = inv(W);

% Return the objective function values after each iteration .
if ret_obj 
    % Define the primal objective function .
    obj = @(X)-log (det (X)) + trace (X*S) + r*sum(abs (X),'all');
    dobj = @(W) log (det(W)) + p ;
    % Initial point is guaranteed to be postive definite (i.e.feasible )
    fvals(1) = obj (X);
    dvals(1) = dobj (W); 
    fprev = fvals(1);
else
    % Track the precision matrix.
    Xprev = X;
end


% Set the initial difference to ensure at least 1 iteration.
rel_diff = 2*tol + 1;

iter = 0; % Track the iteration.
mle_value = [];
while iter < maxiter && rel_diff> tol
% Get the coordinate to update .
    i = mod(iter , p) + 1;
    
    % Get the indices for the block submatrix .
    idx = cat(2,1:i-1,i+1:p);
    
    X11 = X(idx, idx);
    x12 = X(idx, i);
    s12 = S(idx, i);

    % Solve the QP with box constraints. Do a maximum of 100 cycles
    % Initial point is previous optimal value.
    cons = r*ones(p-1, 1);
    dy = QP_box( X11, X11*s12, -cons, cons, x12, 100*p, 1e-6);
    
    % Update the covariance matrix
    w12 = s12 + dy;
    W(idx,1) = w12;
    W(i, idx) =w12'; 


    % Update the precision matrix.
    % Use thresholding to ensure that the update will be sparse.
    x12 = -X11*w12/W(i , i).*(r -abs(dy) < 1e-6);
    x22 = (1 - x12'*w12)/W(i,i);
    X(idx,i) = x12;
    X(i,idx) = x12';
    X(i,i) = x22;
    


    iter = iter+1;
    k_O = r;
    mle_ = mle;
    mle_now =mle_.g(Ydata,Xdata,mu,A,B,X,k_O,k_B);
    mle_value = [mle_value,mle_now];

    % Compute the relative difference for the stopping criteria after
    % every cycle over all coordinates.
    if ret_obj
        % Use relative change in objective function value.
        fvals (iter + 1) = obj (X);
        dvals (iter + 1) = dobj (W);
        if i==p 
            rel_diff = ( fprev - fvals (iter+1))/abs(fprev); 
            fprev = fvals (iter + 1); 
        end
    else
        % Use relative change in precision matrix.
        if i==p
            rel_diff = norm(X - Xprev, 'Fro')/norm(Xprev , 'Fro');
            Xprev = X;
        end
    end
end

X = full(X);
end