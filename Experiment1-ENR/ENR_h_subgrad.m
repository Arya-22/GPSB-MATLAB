function subgrad = ENR_h_subgrad(omega,params)
    % 计算f(ω)的次梯度
    % f(ω) = η2‖ω‖2
    eta2 = params.eta2;  
    
    % 计算二范数 ||x||_2 的次梯度
    % 输入: x - n维向量
    % 输出: g - 次梯度向量
    
    n = length(omega);
    norm_omega = norm(omega, 2);  % 计算二范数
    
    if norm_omega > 1e-6  % x ≠ 0，使用唯一次梯度
        g = omega / norm_omega;
    else
        % x = 0，从单位球中随机选择一个次梯度
        % 方法1: 随机生成单位球内的向量
        u = randn(n, 1);
        g = u / norm(u, 2);
        
        % 方法2: 生成方向随机，长度<=1的向量
        % g = randn(n, 1);
        % g = g / norm(g, 2) * rand;  % 随机长度 [0,1]
    end
    subgrad = g * eta2;
end
