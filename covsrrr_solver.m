function solver = covsrrr_solver
    solver.solveAB = @solveAB;% SRRR
    solver.solveAB_cf = @solveAB_cf;% SRRR
    solver.solveO =@solveO;
end

function [A_best,B_best,mle_in_AB] = solveAB(A,B,Y,X,lr,O,k_1,k_2,max_iter,max_best)

Y_tilde = Y*sqrtm(O);
dim_X = size(X,2);


mle_in_AB = [];
itera = 1; iterb = 1;
obj_val =1e10; obj_best = 1e10; obj_diff = 10000;

while(abs(obj_diff)>0.001 && itera<max_iter && iterb < max_best)
    % 1. Update A for fixed B
    [U,S,V] = svd(Y_tilde'*X*B,"econ");%(q,n)(n,p)(p,r)
    A = U*V';

    cvx_begin quiet
        variable B_cvx(dim_X,lr)
        Reg_B_cvx = 0;
        for i = 1:dim_X
            Reg_B_cvx = Reg_B_cvx + norm(B_cvx(i,:));     
        end 
        minimize(.5*sum(sum_square_abs(transpose(A)*(Y_tilde')-transpose(B_cvx)*(X'))) + k_2*Reg_B_cvx);
        % square_pos(norm(Y,'fro')), or sum_square(vec(Y)), or sum_square_abs(vec(Y)) if Y is complex.
    cvx_end
    B = B_cvx;

    % check whether objective function has converged
    value_now = mle(Y',X',A,B,O,k_1,k_2);
    obj_newval = norm(Y_tilde-X*B*A','fro');
    obj_diff = abs(obj_val - obj_newval);
    obj_val = obj_newval;
    if itera > 1 && obj_val < obj_best
        obj_best = obj_val;
        iter_best = itera;
        B_best = B; A_best = A;
        iterb = 1;
    else
        iterb = iterb+1;
    end
    itera = itera+1;
    
    
end
mle_in_AB=[mle_in_AB,value_now];


fprintf('loop in AB: %6.0f\n',itera);
end

function [A_best,B_best,mle_in_AB] = solveAB_cf(Y,X,lr,O,k_1,k_2,max_iter,max_best)
% if nargin < 5, max_iter = 1e3; end
% if nargin < 6, max_best = 5; end
Y_tilde = Y*sqrtm(O);

dim_X = size(X,2);
B = randn(dim_X,lr);
mle_in_AB = [];

itera = 0; iterb = 0;
obj_val =1e10; obj_best = 1e10; obj_diff = 10000;

while(abs(obj_diff)>0.001 && itera<max_iter && iterb < max_best)
    % 1. Update A for fixed B
    [U,S,V] = svd(Y_tilde'*X*B,"econ");%(q,n)(n,p)(p,r)
    A = U*V';


    for i=1:size(B,1)
        x_i = X(:,i);
        x_k_b_k = 0;
        for k = 1:size(B,1)
            x_k_b_k = x_k_b_k + X(:,k)*B(k,:);
        end
        x_k_b_k = x_k_b_k - X(:,i)*B(i,:);
        r_i = Y*A-x_k_b_k;
        if  1-k_2/(2*norm(x_i'*r_i))<=0
            B(i,:)=0;
        else
            part1 = 1/((x_i')*x_i);
            part2 = (1-k_2/(2*norm(x_i'*r_i)));
            B(i,:)= part1*part2*(x_i')*r_i;
        end      
    end


    

    % check whether objective function has converged
    obj_newval = norm(Y_tilde-X*B*A','fro');
    obj_diff = (obj_val - obj_newval);
    obj_val = obj_newval;
    if itera > 1 && obj_val < obj_best
        obj_best = obj_val;
        iter_best = itera;
        B_best = B; A_best = A;
        iterb = 1;
    else
        iterb = iterb+1;
    end
    itera = itera+1;
    value_now = mle(Y',X',A,B,O,k_1,k_2);
    mle_in_AB=[mle_in_AB,value_now];
end

fprintf('loop in AB: %6.0f\n',itera);
end


function [Omega, mle_in_O] = solveO(O,method,Y,X,A,B,k_1,k_2)
n_sample = size(Y,1);
S = (Y-X*B*(A'))'*(Y-X*B*(A'))/n_sample;
if strcmp('dp',method)% dpglasso(dataX,dataY,A,B,k_1,k_2, S, maxiter, tol)
    [Omega_pre,W,iter_Omega,mle_in_O] = dpglasso(O,X',Y',A,B,k_1,k_2, S, 1e3, 1e-6);
    Omega = full(Omega_pre);
else
   [Omega, W,mle_in_O,cputimeStart] = glasso_1(S,O,k_1,1e-3,Y,X,A,B,k_2);

end
end



