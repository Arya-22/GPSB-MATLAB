%定义次梯度
function grad = ENR_f_subgrad(omega,params)
    % 计算h(ω)的梯度 (h是可微函数)
    % h(ω) = 0.5‖y-Xω‖₂²+eta_1 ||x||_1
    Sigmma_h = params.Sigmma_h;
    Ab = params.Ab; 
    eta1 = params.eta1;
    n = length(omega);
    
    % L1范数的次梯度
    l1_subgrad = zeros(n, 1);
    for i = 1:n
        if omega(i) > 0
            l1_subgrad(i) = 1;
        elseif omega(i) < 0
            l1_subgrad(i) = -1;
        else
            % 在零点，次梯度是[-1, 1]中的任意值
            l1_subgrad(i) = 0; % 可以取0或[-1,1]中的其他值
        end
    end 

    grad =  Sigmma_h * omega - Ab +  l1_subgrad * eta1; 
    
end

