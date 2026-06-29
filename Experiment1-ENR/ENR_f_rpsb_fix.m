function omega_prox = ENR_f_rpsb_fix(v, k, lambda_prox, ki, f_subgrad, params)
    % 使用 FISTA + 回溯线搜索求解临近算子
    % 输入输出不变，内部修正回溯条件与循环保护

    % ========== 1. 提取参数 ==========
    Ab = params.Ab;
    Sigmma_h = params.Sigmma_h;
    x0 = params.x0;
    tau = params.tau;
    theta = params.theta;
    eta1 = params.eta1;
    gpsb = params.gpsb;
    caseki = params.caseki;

    % ========== 2. 确定中心 ==========
    if caseki == 2
        if ki == 2
            center = v;
        else
            center = x0;
        end
    else
        if gpsb == 1
            center = x0;
        else
            center = v;
        end
    end

    % ========== 3. 计算 para ==========
    para = (k - tau(ki) * k^theta(ki)) * Sigmma_h;

    % ========== 4. 构建目标函数的组成部分 ==========
    k_sq = k;
    Ab_vec = k_sq * Ab;
    lambda_l1 = k_sq * eta1;
    inv_lambda = 1 / lambda_prox;

    % 光滑部分 f(x) = 0.5*(x-center)'*para*(x-center) - k_sq*Ab'*x + 0.5*inv_lambda*norm(x-v)^2
    f_smooth = @(x) 0.5 * (x - center)' * para * (x - center) ...
                    - k_sq * Ab' * x + 0.5 * inv_lambda * norm(x - v)^2;
    grad_f = @(x) para * (x - center) - Ab_vec + inv_lambda * (x - v);

    % 完整目标函数（含 l1 范数，用于最终输出）
    fun_val = @(x) f_smooth(x) + lambda_l1 * norm(x, 1);

    % ========== 5. 估计 Lipschitz 常数 ==========
    n = length(v);
    H = para + inv_lambda * speye(n);
    if n <= 1000
        L = max(eig(full(H)));
    else
        L = eigs(H, 1, 'largestabs');
    end
    if L < eps
        L = 1;
    end

    % ========== 6. FISTA 迭代（带回溯线搜索）==========
    max_iter = 500;
    max_backtrack = 50;          % 回溯最大次数
    tol = 1e-4;
    x = v;
    y = x;
    t = 1;
    beta = 0.8;                  % 回溯衰减因子（增大 L）
    L_cur = L;                   % 当前步长倒数（步长 = 1/L_cur）

    for iter = 1:max_iter
        x_old = x;
        grad = grad_f(y);

        % --- 回溯线搜索确定 L_cur ---
        for bt = 1:max_backtrack
            % 梯度步 + 软阈值（近端步）
            z = y - (1 / L_cur) * grad;
            x_new = sign(z) .* max(abs(z) - lambda_l1 / L_cur, 0);

            % 检查光滑部分的二次上界
            rhs = f_smooth(y) + grad' * (x_new - y) + (L_cur / 2) * norm(x_new - y)^2;
            if f_smooth(x_new) <= rhs + 1e-12
                break;   % 接受当前步长
            end
            % 不满足，增大 L（减小步长）
            L_cur = L_cur / beta;
        end
        if bt == max_backtrack
            warning('Backtracking reached max iterations, L_cur = %e', L_cur);
        end

        x = x_new;
        L = L_cur;   % 更新 L 用于下次迭代

        % --- FISTA 加速 ---
        t_old = t;
        t = (1 + sqrt(1 + 4 * t^2)) / 2;
        y = x + (t_old - 1) / t * (x - x_old);

        % --- 收敛判断 ---
        rel_change = norm(x - x_old) / (norm(x_old) + eps);
        if rel_change < tol
            break;
        end
    end

    if iter == max_iter
        warning('FISTA reached max iterations (%d) without converging. rel_change = %e', max_iter, rel_change);
    end

    omega_prox = x;
end