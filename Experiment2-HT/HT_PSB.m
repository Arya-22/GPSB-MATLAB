
function [k, Fv,gapFF,gapx,x_aver,cpu_time,F_value] = HT_PSB(F,f_subgrad, h_subgrad, f_prox, params)
    
    % 初始化
    x0 = params.x0;            % 初始点
    max_iter = params.max_iter;
    tol = params.tol;         % 收敛容忍度
    x_true =  params.original_signal;% 收敛容忍度
    ki = params.ki;
    % A = params.A;
    % n = params.n;
    % m = params.m;
    % b = params.b; 
    gamma = params.gamma(ki);
    beta_k = 1;
    s_k=0;
    x = x0;
    m=length(x);
    history = zeros(max_iter,m);
    Fv = zeros(max_iter,1);
    gapFF = zeros(max_iter,1);
    gapx = zeros(max_iter,1);
    cpu_time = zeros(max_iter,1);
    theta_0 = params.theta_0; 
    % Sigmma_h = params.Sigmma_h;
    Ab = params.Ab;

% 预分配表格
fprintf('%-12s %-10s %-12s %-12s %-12s %-10s\n', 'ki', '迭代次数','参数Fv', 'gapFF', '结果m', 'CPU时间');

    % 主迭代循环
    for k = 1:max_iter
        tic 
        % 1. 计算次梯度
        subgrad = h_subgrad(x); 
        stepsize = 1; % 初始步长系数
        lambda_prox = 1/(gamma*beta_k);
        
        if ki==1 %MPSB1
                s_k = s_k + stepsize*subgrad;
                v = x0 - s_k/(gamma*beta_k); 
                x_next = f_prox(v,k,lambda_prox,ki); 
%             s_k = s_k + stepsize*subgrad;
%             x_next = x0 - k * lambda_prox * Ab/n - s_k * lambda_prox;
            elseif ki==2 %MPSB2 RPSB 初始参数 \tau theta 不同
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - s_k/(gamma*beta_k);
                v = z_k;
                x_next = f_prox(v,k,lambda_prox,ki);  
%             x_next = - k * lambda_prox * Ab/n + (eye(m,m)- k*lambda_prox*Sigmma_h)*(x0  - s_k * lambda_prox);
            elseif ki==3 %MPSB2 RPSB 初始参数 \tau theta 不同
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - s_k/(gamma*beta_k);
                v = z_k;
                x_next = f_prox(v,k,lambda_prox,ki); 
            elseif ki == 4 %DA
                s_k = s_k + stepsize*(subgrad + f_subgrad(x));
                x_next = x0-s_k/(gamma*k^0.5);
            elseif ki==5 %PSM
                v = x - subgrad*beta_k/(gamma*k); 
                x_next = f_prox(v,k,lambda_prox,ki);
            elseif ki==6 %SM 
                x_next = x - (subgrad + f_subgrad(x))/(gamma*k^0.5);
        end
            gapx0 = 2 * norm(x0-x_true,2)^2;
                diff = x - x0;
                norm_diff = norm(diff, 2);
        
        %3.如果存在约束集，要加一步投影        
                % 计算半径
                radius = sqrt(2 * gapx0);
                
                % 判断是否在集合内
                if norm_diff <= radius
                    % 如果已经在集合内，保持不变
                    x = x;
                else
                    % 投影到球面上
                    x = x0 + (radius / norm_diff) * diff;
                end
        

        %4. 记录与收敛检查
        history(k, :) = x;  
        x_aver =history(1:k,:)' * ones(k,1)/k;
        cpu_time(k) = toc;
        if k>1
            Fv(k, :) = F(x_aver);
            x_aver_old =history(1:(k-1),:)' * ones(k-1,1)/(k-1);
            gapFF(k,:) = abs(F(x_aver)-F(x_aver_old)); 
            if gapFF(k,:) < gapFF(k-1,:)
                gapFF(k,:)  = gapFF(k,:);
            else
                gapFF(k,:) = gapFF(k-1,:);
            end
            gapx(k, :)  = norm(x-x0,2); 
            % if k > max_iter || 0.5 * gapx(k,:)^2 > gapx0
            if gapFF(k,:) < tol 
            % if F(x_aver) < F(x_true) && 0.5 * gapx(k,:)^2 > gapx0
            % if sum(cpu_time(1:k)) > m/30
                fprintf('%-10d %-10d %-10d %-10d %-12.4f %-10.4f\n', ki, k,Fv(k), gapFF(k,:), m, cpu_time(k));
                break;
            end
        else 
            gapFF(1,:) = abs(F(x0)-F(x_true))/F(x_true);  
        end 

      
        % 5. 更新迭代点
        x = x_next; 
            beta_k= k^theta_0(ki); 
    fprintf('%-10d %-10d %-10d %-12.4f%-12.4f %-10.4f\n', ki, k, Fv(k), gapFF(k,:), m, cpu_time(k));
    end
    F_value = F(x_aver);
    % errorF  = sum(A*x_aver-b)/n;
    % res = 1- norm((A*x_aver-b))/norm(b-mean(b));
end