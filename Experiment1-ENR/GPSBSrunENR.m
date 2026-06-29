    clear all
    clc
    % colors =  [...
    %     'r',       ... % 红色 (Red)
    %     'g',       ... % 绿色 (Green)
    %     'b',       ... % 蓝色 (Blue)
    %     'm',       ... % 洋红 (Magenta)
    %     'c',       ... % 青色 (Cyan)
    %     'k'
    % ];
    % linestyles = {...
    %     '-',       ... % 实线 (Solid)
    %     '--',      ... % 虚线 (Dashed)
    %     ':',       ... % 点线 (Dotted)
    %     '-.',      ... % 点划线 (Dash-Dot)
    %     '--',      ... % 虚线 (Dashed)
    %     '-',       ... % 实线 (Solid)
    %     ':',       ... % 点线 (Dotted)
    %     '-.',      ... % 点划线 (Dash-Dot)
    %     '--',      ... % 虚线 (Dashed)
    %     '-'         ... % 实线 (Solid)
    % };
    % markers = {...
    %     'o',       ... % 圆圈 (Circle)
    %     '+',       ... % 加号 (Plus)
    %     '*',       ... % 星号 (Asterisk)
    %     's',       ... % 方块 (Square)
    %     'd',       ... % 菱形 (Diamond)
    %     '^',       ... % 上三角 (Up Triangle)
    %     'v',       ... % 下三角 (Down Triangle)
    %     '>',       ... % 右三角 (Right Triangle)
    %     '<',       ... % 左三角 (Left Triangle)
    %     'p'         ... % 五角星 (Pentagram)
    % };

    r=1;%figure 排序
for p = 1% p=1,2,3,4的选取是方便caseki=3时参数的循环，对其他情况单选不影响
   for caseki = 3 
    for gpsb = 1:2
    params.gpsb = gpsb;%参数对比时运行gpsb1选1，运行gpsb2选2，算法对比时，不影响不用管
    if caseki == 0 %tau参数灵敏度对比
        tau = [0,0.3,0.5,0.7,0.9];
        k_num = length(tau);
        theta = 0.5*ones(1,k_num);
        theta_0 = 0.9*ones(1,k_num);
        num = 1:k_num; 
        gamma = 5e1*ones(1,length(tau));
%       gamma = [2e-1,7e-1,5e-1,2e-1,5e0];
    elseif caseki == 1 %theta参数灵敏度对比
        theta = [0,0.1,0.3,0.5,0.8];
        k_num = length(theta);
        num = 1:length(theta);%运行theta的参数/对比RPSB,GPSB1,IPSB,GPSB2
        tau = 0.2*ones(1,k_num);
        theta_0 = 0.9*ones(1,k_num);
      gamma = 5e1*ones(1,length(tau));
        % gamma = [7e1,2e1,5e-1,2e-1,5e-1];
    elseif caseki == 2  %算法对比
       tau = [0.1,0.1,0,1,0,2];% GPSB1,GPSB2,RPSB,DA,PSM,SM
       theta = [0.3,0.4,1,0.5,1,2];%t_value_1234  = tau.*k.^(theta);
       num = 1:length(tau);
       theta_0 = 0.4*ones(1,length(tau));
        gamma = [5e1,5e1,5e1,5e3,5e3,5e3];
      % gamma = 5e1*ones(1,length(tau));
    elseif caseki == 3 %theta0参数灵敏度对比 
        if p == 1
            theta_0 = [0.1,0.3,0.5,0.7,0.9];
            k_num = length(theta_0); 
            tau = 0.1*ones(1,k_num);
            theta = 0.6*ones(1,k_num);  
            gamma = 1e2*ones(1,length(tau));
        elseif p==2 
            theta_0 = [0.1,0.3,0.5,0.7,0.9];
            k_num = length(theta_0); 
            tau = 0.6*ones(1,k_num);
            theta = 0.1*ones(1,k_num);   
            gamma = 1e2*ones(1,length(tau));%GPSB2
        elseif p==3
            theta_0 = [0.1,0.3,0.5,0.7,0.9];
            k_num = length(theta_0); 
            tau = 0.5*ones(1,k_num);
            theta = 0.5*ones(1,k_num);  
            gamma = 1e2*ones(1,length(tau));
        elseif p==4
            theta_0 = [0.1,0.3,0.5,0.7,0.9];
            k_num = length(theta_0); 
            tau = 0.0*ones(1,k_num);
            theta = 0.0*ones(1,k_num);  
            gamma = 1e2*ones(1,length(tau));
        end
        num = 1:length(theta_0);
    else
        fprintf('1是参数theta对比，2是固定参数下的算法对比,3是参数theta_0对比');
    end
    params.caseki = caseki;
    data = data2_Friedman();
    A = data.A;
    b = data.b;
    original_signal = data.x;
    m = data.m;
    n = data.n;  
    H=data.H; 
    SNR = data.SNR; 
    % 生成数据    
    rng('default');
    rng(1);
    % 预分配结果数组
    max_iter = 300;
    eta1 = 0.01;
    eta2 = 0.02;
    k_num = length(num); 
    gapx_rpsb = zeros(max_iter,length(m),k_num); 
    gapF_rpsb = zeros(max_iter,length(m),k_num); 
    cpu_time_rpsb = zeros(max_iter,length(m),k_num);
    k_rpsb =zeros(length(m),k_num); 
    errorF =zeros(length(m),k_num); 
    res =zeros(length(m),k_num);  
    cpu_time_sum = zeros(1,k_num);
    
    
    params.max_iter = max_iter ;
    params.tol = 1e-3; 
    params.dim = m; 
    params.num = k_num;
    params.tau = tau;
    params.theta = theta;
    params.A = A;
    params.b = b;
    params.n = n;
    params.eta1 = eta1;
    params.eta2 = eta2;
    x_aver =zeros(m,length(m),k_num);
    for i = 1:length(m)         
            % 定义函数句柄
            
            rng(1);
            params.x0 =zeros(m,1);
            params.Ab = A'*b;       
            Sigmma_h = A' * A;
            params.Sigmma_h = Sigmma_h; 
            f_subgrad = @(x) ENR_f_subgrad(x,params);
            h_prox = @(v,lam) ENR_h_prox(v,lam); 
            F = @(x) (1/2) * norm(A * x - b,2)^2 + eta1*norm(x,1) + eta2*norm(x,2);  
            h_subgrad = @(x) ENR_h_subgrad(x, params);         
            params.original_signal = original_signal;

            
            f_rpsb = @(v, k,lambda_prox,ki) ...
                ENR_f_rpsb_fix(v, k,lambda_prox,ki,f_subgrad,params);

        for ki = 1:k_num
            params.gamma = gamma;
            params.theta_0 = theta_0(ki);%区间为(0,1)
            params.ki = ki;
            [k_rpsb(i,ki),gapF_rpsb(:,i,ki),gapx_rpsb(:,i,ki),x_aver(:,i,ki), cpu_time_rpsb(:,i,ki),errorF(i,ki),res(i,ki),x_old] = ...
                PSB(F,f_subgrad,h_subgrad,f_rpsb,params);
                    cpu_time_sum(i,ki) = sum(cpu_time_rpsb(:,i,ki)); 
        end
    end
 

    % % 调数据 单独画图 ERROE epsilon_k
    % % k_num = length(tau);
    % colors = lines(k_num); % 获取k_num种不同的颜色
    % linestyles = { '-.','--',':','-',':','--'}; % 定义线型
    % markers = {...
    %     'o',       ... % 圆圈 (Circle)
    %     '+',       ... % 加号 (Plus)
    %     '>',       ... % 星号 (Asterisk)
    %     'square',       ... % 方块 (Square)
    %     '<',       ... 
    %     '*',       ... % 星号 (Asterisk)
    %     };
    % 
    % 
    % for ii = 1:length(m)
    %     % 创建并设置gap图形 
    %     figGap = figure(r);
    %     r = r+1;
    %     clf(figGap);
    %     set(figGap,'Name',['m=',num2str(m(ii)),' n=',num2str(n(ii))]);
    %     hold on;
    %     % legend_label = cell(1,k_num);
    %     for ki = 1:k_num % 假设有k_num种画法
    %         plot(2:5:k_rpsb(ii,ki), log(gapF_rpsb(2:5:k_rpsb(ii,ki),ii,ki)),...
    %             linestyles{ki},'Color',colors(ki,:),...
    %             'Marker',markers{ki},...
    %             'MarkerSize',8,...
    %             'LineWidth',2.5,...
    %             'MarkerIndices', 2:5:length(2:5:k_rpsb(ii,ki))); 
    %         % 根据caseki生成图例标签
    %         if caseki == 2 
    %             method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
    %             legend_label  = method_names; 
    % %             xlim([3,90]); 
    % %             ylim([0,3]); 
    %         elseif caseki == 0
    %             legend_label{ki} = sprintf("$\\tau= $%.2f", tau(ki));  
    %             title(['$(\theta_0,\theta) =$','(',num2str(theta_0(ki)),',',num2str(theta(ki)),')'], 'Interpreter','latex', 'FontSize', 12); 
    %         elseif caseki == 1
    %             legend_label{ki} = sprintf("$\\theta= $%.2f", theta(ki));  
    %             title(['$(\theta_0,\tau) = $','(',num2str(theta_0(ki)),',',num2str(tau(ki)),')'], 'Interpreter','latex', 'FontSize', 12);
    % %             xlim([3,60]); 
    % %             ylim([0,400]); 
    %         elseif caseki == 3
    %             legend_label{ki} = sprintf("$\\theta_0= $%.2f", theta_0(ki));
    %             title(['$(\tau,\theta) = $','(',num2str(tau(ki)),',',num2str(theta(ki)),')'], 'Interpreter','latex', 'FontSize', 12);
    % %             ylim([0,230]); 
    %         else
    %            fprintf('caseki不等于1，2，3')
    %         end
    %     end
    %     hold off;
    % 
    %     % 循环结束后添加图例
    %     legend(legend_label, 'Interpreter', 'latex', 'FontSize', 18);
    %     xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 16);
    %     ylabel('Error Tolerance $log (\epsilon_k) $','Interpreter','latex', 'FontSize', 16);
    %     grid on;
    %     set(gca, 'FontSize', 16, 'GridAlpha', 0.3);
    %     set(gcf, 'Position', [100, 100, 800, 600]);
    % end

%% 
% 调数据 单独画图 MSE
k_num = length(tau);
colors = lines(k_num); % 获取k_num种不同的颜色
linestyles = { '-.','--',':','-',':','--'}; % 定义线型
markers = {...
    'o',       ... % 圆圈 (Circle)
    '+',       ... % 加号 (Plus)
    '>',       ... % 星号 (Asterisk)
    's',       ... % 方块 (Square)
    '<',       ... 
    '*',       ... % 星号 (Asterisk)
    }; 
% for k = 1:max(size(cpu_time_rpsb)) %以时间为下坐标
%     for ii = 1:length(m)
%         for ki = 1:k_num % 假设有k_num种画法
%             cpu_time_rpsb_e(k,ii,ki)= sum(cpu_time_rpsb(1:k,ii,ki));
%         end
%     end
% end
for ii = 1:length(m)
    % 创建并设置gap图形 
    figGap = figure(r); 
    clf(figGap);
    set(figGap,'Name',['EnrMse','m=',num2str(m(ii)),' n=',num2str(n(ii))]);
    hold on;
    % legend_label = cell(1,k_num);
    
    for ki = 1:k_num % 假设有k_num种画法 
%         plot(cpu_time_rpsb_e(2:5:k_rpsb(ii,ki),ii,ki),log(gapx_rpsb(2:5:k_rpsb(ii,ki),ii,ki)),...
%             linestyles{ki},'Color',colors(ki,:),...
%             'Marker',markers{ki},...
%             'MarkerSize',8,...
%             'LineWidth',2.5,...
%             'MarkerIndices',2:5:length(2:5:k_rpsb(ii,ki))); 
        plot(2:5:k_rpsb(ii,ki),log(gapx_rpsb(2:5:k_rpsb(ii,ki),ii,ki)),...
            linestyles{ki},'Color',colors(ki,:),...
            'Marker',markers{ki},...
            'MarkerSize',8,...
            'LineWidth',2.5,...
            'MarkerIndices',2:5:length(2:5:k_rpsb(ii,ki))); 
        % 根据caseki生成图例标签
        if caseki == 2 
            method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
            legend_label  = method_names; 
%             xlim([3,90]); 
%             ylim([0,3]); 
        elseif caseki == 0
            legend_label{ki} = sprintf("$\\tau= $%.2f", tau(ki));  
            title(['$(\theta_0,\theta) =$','(',num2str(theta_0(ki)),',',num2str(theta(ki)),')'], 'Interpreter','latex', 'FontSize', 12); 
        elseif caseki == 1
            legend_label{ki} = sprintf("$\\theta= $%.2f", theta(ki));  
            title(['$(\theta_0,\tau) = $','(',num2str(theta_0(ki)),',',num2str(tau(ki)),')'], 'Interpreter','latex', 'FontSize', 12);
%             xlim([3,60]); 
%             ylim([0,400]); 
        elseif caseki == 3
            legend_label{ki} = sprintf("$\\theta_0= $%.2f", theta_0(ki));
            title(['$(\tau,\theta,\gamma) = $','(',num2str(tau(ki)),',',num2str(theta(ki)),',',num2str(gamma(ki)),')'], 'Interpreter','latex', 'FontSize', 12);
%             ylim([0,230]); 
        else
           fprintf('caseki不等于1，2，3')
        end
    end
    hold off;

    % 循环结束后添加图例
    legend(legend_label, 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 16);
    ylabel('log MSE','Interpreter','latex', 'FontSize', 16);
    grid on;
    set(gca, 'FontSize', 16, 'GridAlpha', 0.3);
    set(gcf, 'Position', [800, 200, 800, 600]);
end
%%
% 存图   
        saveDir = 'D:\20260626'; 
            % 确保文件夹存在，如果不存在则创建
            if ~exist(saveDir,'dir')
                mkdir(saveDir); 
            end  
        % 获取所有图形
        figures = findall(groot,'Type','figure');
        for idx = 1:length(figures)
        % 切换到当前图形
        figure(figures(idx));    
        % 获取标题（尝试从坐标轴获取）
            if caseki == 1
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptheta'];
            elseif caseki == 2
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cpms'];
            elseif caseki == 3
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptheta0','_p',num2str(p)];
            elseif caseki == 0
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptau'];
            else
                fprintf('WRONG NAME');
            end
        % 处理非法字符
        fileName = regexprep(fileName,'[\\/:*?"<>|]','_');

        fileName = [fileName,['_gpsb',num2str(gpsb),'_caseki',num2str(caseki)]];
        % 完整文件路径
        fullPath = fullfile(saveDir,fileName);   
        %     保存为EPS格式（矢量图）
        print(figures(idx),[fullPath '.eps'],'-depsc','-r600');   
        % 保存为FIG格式
        savefig(figures(idx),[fullPath '.fig']);
        end


   
        % % 存数据 
        saveDir = 'D:\20260626';
     % 确保文件夹存在，如果不存在则创建
            if ~exist(saveDir,'dir')
                mkdir(saveDir); 
            end  
        for ii = 1:length(m)
            % 拼接完整文件路径并保存        
            if caseki == 1
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptheta'];
            elseif caseki == 2
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cpms'];
            elseif caseki == 3
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptheta0','_p',num2str(p)];
            elseif caseki == 0
                fileName = ['m',num2str(m(ii)),'n',num2str(n(ii)),'_Cptau'];
            else
                fprintf('WRONG NAME');
            end
            % 处理非法字符
        fileName = regexprep(fileName,'[\\/:*?"<>|]','_');
        fileName = [fileName,['gpsb',num2str(gpsb)]]; 

        fullFileName = fullfile(saveDir, fileName); % fullfile 能正确处理不同系统的路径分隔符

            save(fullFileName,"A","b","cpu_time_sum","caseki","eta1","eta2","errorF","gamma","fileName",...
                "gapF_rpsb","gapx_rpsb","k_rpsb","m","n", "tau","theta","theta_0","x_aver","gapF_rpsb",...
                "original_signal","res","colors","legend_label","linestyles","markers","gpsb");  % 保存所有工作区变量
        end 
        r=r+1;
    end 
   end
end
 


