function x_prox = HT_f_sigmmaf_prox(v, k,lambda_prox,ki,params)
    % $f(x)=\gamma\|x\|_1+\frac{1}{n} {\|Ax-b\|_1}; \\h(x)=-\frac{1}{2n}  \|Ax - b\|^2_2$
    %计算h(ω)的临近算子
    A = params.A;
    b = params.b;
    Sigmma_h = params.Sigmma_h;
    Ab = params.Ab;
    tau = params.tau;
    theta = params.theta;
    [n,~] = size(A); 
    x0 = params.x0;
    Ab = params.Ab;
    nb = params.nb;
    para = (tau(ki)*k^theta(ki)) * Sigmma_h;
    nx0 = x0'/(2*n)* para * x0; 
 
if ki == 3 %RPSB初始点v不同
    fun = @(x) k/2 * (x'*Sigmma_h * x) + k*nb  + k/n * x'* Ab   + 0.5 * (1/lambda_prox) * norm(x - v)^2;
    grad_fun = @(x)  k/n * A' * b  +  (1/lambda_prox) * (x - v) ;
elseif ki == 5 %PSM 初始点v不同
    fun = @(x)  k/2 * (x'*Sigmma_h * x) + k*nb  + k/n * x'* Ab   + 0.5 * (1/lambda_prox) * norm(x - v)^2;
    grad_fun = @(x)  k/n * A' * b  +  (1/lambda_prox) * (x - v) ;
elseif ki == 1 %GPSB1  初始点不同 
    fun = @(x) k*nb + k/n * x'* Ab + nx0 +  x'/(2*n)* para * x0 + 0.5 * (1/lambda_prox) * norm(x - v)^2;
    grad_fun = @(x) k/n * Ab + 2 * k * para * v +  (1/lambda_prox) * (x - v) ; 
elseif ki==2 %GPSB2
    fun = @(x) k*nb + k/n * x'* Ab -  (v'/(2*n)-  x'/(2*n))* para * v + 0.5 * (1/lambda_prox) * norm(x - v)^2;
    grad_fun = @(x) k/n * Ab + 2 * k * para * v +  (1/lambda_prox) * (x - v) ;
end
 

% MINFUNCTION 初始点（从 v 开始）
x0 = v;

    options = [];
    options.Method = 'lbfgs';
    options.MaxIter = 1000;
    options.Display = 'off';
    options.TolFun = 1e-8; 
    x_opt = minFunc(@objective, x0, options);
    
    function [f, g] = objective(x)
        f = fun(x);
        g = grad_fun(x);
    end


% 输出结果
x_prox = x_opt;

end

