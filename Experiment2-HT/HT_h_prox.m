%定义临近算子
function omega_prox = HT_h_prox(v, params)
    % 计算f(ω)的临近算子
    %$f(x)=\gamma\|x\|_1+\frac{1}{n} {\|Ax-b\|_1}; \\h(x)=-\frac{1}{2n}  \|Ax - b\|^2_2$
     
    A = params.A;
    b = params.b;
    gamma = params.gamm;
    % 复合L1范数函数的临近算子
    % f(x) = gamma * ||x||_1 + (1/n) * ||A*x - b||_1
    
    rho = 1.0; 
    max_iter = 1000; 
    tol = 1e-6; 
    
    [n, m] = size(A);
    y = zeros(n, 1);
    u = zeros(n, 1);
    
    % 预计算用于x子问题的矩阵
    I = eye(m);
    if m <= 1000  % 小规模问题直接求逆
        M = I + rho * (A' * A);
        solver = @(rhs) M \ rhs;
    else  % 大规模问题使用共轭梯度
        solver = @(rhs) pcg(@(x) x + rho * A' * (A * x), rhs, 1e-8, 100);
    end
    
    for iter = 1:max_iter
        % x子问题
        c = b + y - u;
        rhs = v + rho * A' * c;
        x_new = solver(rhs);
        
        % 软阈值（L1正则化部分）
        x_new = sign(x_new) .* max(abs(x_new) - gamma, 0);
        
        % y子问题
        Axb = A * x_new - b;
        y_new = sign(Axb + u) .* max(abs(Axb + u) - 1/(n*rho), 0);
        
        % 对偶变量更新
        u = u + Axb - y_new;
        
        % 收敛检查
        primal_res = norm(A * x_new - b - y_new);
        dual_res = rho * norm(A' * (y_new - y));
        
        if primal_res < tol && dual_res < tol
            break;
        end
        
%         x = x_new;
        y = y_new;
    end
    
    omega_prox = x_new;
end