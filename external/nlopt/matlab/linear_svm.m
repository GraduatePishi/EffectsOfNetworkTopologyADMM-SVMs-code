function [xave, history] = linear_svm(A, lambda, p, rho, alpha, fileID3)
t_start = tic;
QUIET    = 0;
MAX_ITER = 10000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

[m, n] = size(A);
p=ones(1,length(p)); %testing
N = max(p);
% group samples together
for i = 1:N,
    tmp{i} = A(p==i,:);
end
A = tmp;


x = zeros(n,N);
z = zeros(n,N);
u = zeros(n,N);



if ~QUIET
    fprintf(fileID3,'%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end


for k = 1:MAX_ITER

	% x-update
%x_test = zeros(3,1);
    for i = 1:N,
        %x_var = x(:,max(i-1,1));
        %tmp = sum(pos(A{i}*x_var + 1)) + rho/2*sum_square(x_var - z(:,i) + u(:,i))
        cvx_begin quiet
            variable x_var(n)
             minimize ( sum(pos(A{i}*x_var + 1)) + rho/2*sum_square(x_var - z(:,i) + u(:,i)) )
        cvx_end
        
%     A{i}
%     Ax=A{i}*x_test
%     Ax1=A{i}*x_test + 1
%     posAx1=pos(A{i}*x_test + 1)
%     test1 =  sum(pos(A{i}*x_test + 1)) 
%     test2= rho/2*sum_square(x_test - z(:,i) + u(:,i)) 
%     test= test1+ test2
        x(:,i) = x_var;
    end
  
    xave = mean(x,2);

    % z-update with relaxation
    zold = z;
    x_hat = alpha*x +(1-alpha)*zold;
    z = N*rho/(1/lambda + N*rho)*mean( x_hat + u, 2 );
    z = z*ones(1,N);

    % u-update
    u = u + (x_hat - z);

    % diagnostics, reporting, termination checks
    
    history.objval(k)  = objective(A, lambda, p, x, z);

    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);

    if ~QUIET
        fprintf(fileID3, '%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
        
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end
end

if ~QUIET
    toc(t_start);
end
end

function obj = objective(A, lambda, p, x, z)
    obj = hinge_loss(A,x) + 1/(2*lambda)*sum_square(z(:,1));
end

function val = hinge_loss(A,x)
    val = 0;
    for i = 1:length(A)
        val = val + sum(pos(A{i}*x(:,i) + 1));
    end
end

