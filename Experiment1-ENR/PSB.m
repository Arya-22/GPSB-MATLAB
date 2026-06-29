
function [k, gapF,gapx,x_aver,cpu_time,errorF,res,x_old] = PSB(F,f_subgrad, h_subgrad, f_prox, params)
    
    % 初始化
    x = params.x0;            % 初始点
    max_iter = params.max_iter;
    tol = params.tol;         % 收敛容忍度
    eta1 =  params.eta1;% 收敛容忍度
    ki = params.ki;
    gamma = params.gamma(ki);% 初始步长系数
    caseki = params.caseki;
    A = params.A;
    b = params.b;
    n = params.n;
    Ab = params.Ab; 
    Sigmma_h = params.Sigmma_h; 

    hbeta_k=1;
    s_k=0;
    x0 = x;
    m=length(x);
    history = zeros(max_iter,m);
    gapF = zeros(max_iter,1); 
    gapx = zeros(max_iter,1);
    cpu_time = zeros(max_iter,1); 
    theta_0 = params.theta_0;   
    tau = params.tau;
    theta = params.theta;

% 预分配表格
fprintf('%-12s %-10s %-12s %-12s %-12s %-10s\n', 'ki', '迭代次数','参数gapx', '参数gapF', '结果m', 'CPU时间');

    % 主迭代循环
    for k = 1:max_iter
        tic
        % 1. 计算次梯度
        subgrad = h_subgrad(x); 
        stepsize = 1;% 初始步长系数
        
        lambda_prox = 1/(gamma*hbeta_k);
%         caseki = 0 对比tau，caseki = 1 对比theta, caseki = 2对比算法，caseki=3对比theta0
        if caseki == 2  
            if ki==1%GPSB1 
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - lambda_prox * s_k;
                v=z_k;
                x_next = f_prox(v,k,lambda_prox,ki);   
            elseif ki==2%GPSB2  
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - lambda_prox * s_k;
                v=z_k;
                x_next = f_prox(v,k,lambda_prox,ki);  
            elseif ki==3%RPSB 初始参数 \tau theta 不同
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - s_k *lambda_prox;
                v = z_k;
                x_next = f_prox(v,k,lambda_prox,ki); 
            elseif ki == 4 %DA
                s_k = s_k + stepsize*(subgrad + f_subgrad(x));
                x_next = x0 - lambda_prox * s_k;
            elseif ki==5%PSM
                v = x - subgrad*(1/(gamma*hbeta_k)); 
                x_next = f_prox(v,k,lambda_prox,ki);
            elseif ki==6%SM 
                x_next = x - (subgrad + f_subgrad(x))*lambda_prox;
            end
        else  
            %GPSB1
                s_k = s_k + stepsize*subgrad;
                z_k = x0 - lambda_prox * s_k;
                v=z_k;
                x_next = f_prox(v,k,lambda_prox,ki); 
        end
        
%         
        %3.如果存在约束集，要加一步投影
        

        % 4. 记录与收敛检查
        history(k, :) = x; 
        x_aver =history(1:k,:)' * ones(k,1)/k;
        if k >1
            x_old =history(1:(k-1),:)' * ones((k-1),1)/(k-1); 
              gapx(k, :) = norm(A*x_aver-b)^2/n;%MSE
              gapF(k, :) = norm(A*x_aver-A*x_old)+norm(x_aver-x_old);
            if gapF(k, :)  < tol
                fprintf('%-10d %-10d %-10d %-10d %-12.4f %-10.4f\n', ki, k,gapx(k), gapF(k), m, cpu_time(k));
                break;
            end
        end
%         if k>2
%             if gapF(k, :)+ gapF(k-1, :) < 2*tol 
% %             if sum(cpu_time) > m/5e2
%                 fprintf('%-10d %-10d %-10d %-12.4f %-10.4f\n', ki, k,gapF(k), m, cpu_time(k));
%                 break;
%             end
%             x_old =history(1:(k-1),:)' * ones((k-1),1)/(k-1);
%             gapF(k, :) = norm(A*history(k, :)'-A*history(k-1,:)')+norm(history(k, :)-history(k-1, :));
%             gapx(k,:) = norm(A*x_aver-b)/n;
%         else
%             x_old = x0;
%             gapF(k, :) = norm(A*history(k, :)'-A*x0)+norm(history(k, :)'-x0);
%             gapx(k,:) =  norm(A*x_aver-b)/n; 
%         end
      
        % 5. 更新迭代点
        x = x_next; 
        hbeta_k= k^theta_0; 
        % hbeta_k= k^theta_0-(tau(ki)*k^theta(ki))*sqrt(min(abs(eig(params.Sigmma_h))))/gamma;

       cpu_time(k) = toc;
    fprintf('%-10d %-10d %-12.4f %-12.4f %-12.4f %-10.4f\n', ki, k, gapx(k),gapF(k), m, cpu_time(k));
    end
    errorF  =norm(A*x_aver-b)^2/n;
    res = 1- norm(A*x_aver-b)/sum((b-mean(b)).^2);
end