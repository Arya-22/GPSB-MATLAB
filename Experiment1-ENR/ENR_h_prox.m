%定义临近算子
function omega_prox = ENR_h_prox(v, t, eta1, eta2)
    % 计算f(ω)的临近算子
    % prox_{t*f}(v) = argmin_ω { f(ω) = η2‖ω‖2+ 1/(2t)‖ω-v‖₂²}
    % 
    
    % 解析解: 对于弹性网正则化，可以分解为坐标-wise 的软阈值
    % ω_i = sign(v_i) * max(|v_i| - t*η₁, 0) / (1 + 2t*η₂)
    %     omega_prox = sign(v) .* max(abs(v) - t * eta1, 0) / (1 + 2 * t * eta2);
    
%     omega_prox = sign(v) .* max(abs(v) - t * eta1, 0) / (1 + 2 * t * eta2);
end