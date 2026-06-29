
 
function l1_prox = l1prox(z,lambda)

    % 1范数的临近算子（软阈值算子）
    % 输入:
    %   v - 输入向量/矩阵
    %   lambda - 正则化参数 (λ > 0)
    % 输出:
    %  临近算子结果
    
    if nargin < 2 || isempty(lambda)
        lambda = 1;
    end
    
    if lambda <= 0
        error('lambda must be positive');
    end
    
    % 软阈值操作
    l1_prox = sign(z) .* max(abs(z) - lambda, 0);  
end