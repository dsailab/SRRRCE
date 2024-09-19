function solver = SRRRCE_solver
solver.solve_mu = @solve_mu;
solver.solve_A = @solve_A;
solver.solve_B = @solve_B;
solver.solve_O = @solve_Omega_mm;
end

function mu = solve_mu(Y,X,mu,A,B,O,n)
mu = (1/n)*(Y-A*B'*X)*ones(n,1);
end

function A=solve_A(Y,X,mu,A,B,O,n)

BXXB=(B')*X*(X')*B;
lambda=max(eig(BXXB));

Q_A = A - (1/lambda)*A*B'*X*(X')*B + (1/lambda)*(Y-mu*ones(1,n))*(X')*B;

[U,S,V]=svd(O*Q_A,"econ");%SVD econ 
A=U*V';
end

function B=solve_B(Y,X,mu,A,B,O,n,k_B)
lambda = max(eig(A'*O*A))*max(eig(X*X'));
QB = B - (1/lambda)*(X*X')*B*A'*O*A ...
    +(1/lambda)*X*((Y-mu*ones(1,n))')*O*A;



for i=1:size(B,1)
    if  norm(QB(i,:))<= k_B*(1/lambda)
        B(i,:)=0;
    else
        B(i,:)=(1-k_B*(1/(lambda*norm(QB(i,:)))))*QB(i,:);
    end      
end
end

function O=solve_Omega_mm(Y,X,mu,A,B,O,n,k_O)
S=(1/n)*(Y-mu*ones(1,n)-A*(B)'*X)*(Y-mu*ones(1,n)-A*(B)'*X)';


for j=1:size(S,1)
    % permutation
    idx = 1:size(S,1); idx(j) = [];
    O11 = O(idx, idx);
    o12_last = O(idx, j);
    s22 = S(j,j);
    s12 = S(idx, j);


    
    % surrogate direct estimation
    inv_O11 = inv(O11); % inverse of Omega11
    lambdaO = max(eig(inv_O11)); % max eigenvalue of the inverse of Omega11

    neg_qwi_t = (1/s22)*(1/lambdaO)*(s22*(o12_last')*(inv_O11-lambdaO*eye(size(S,1)-1))+s12');
    qwi =  -(neg_qwi_t');


    o12_est = ell1_prox(qwi,k_O/(s22*lambdaO));

    o22_est=(1/s22)+(o12_est')*inv_O11*o12_est;

    O(idx,j) = o12_est;
    O(j,idx) = o12_est';
    O(j,j) = o22_est;
end
end


