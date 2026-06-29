%%
clearvars
clc

% 调数据 单独画图
% 预先定义图形参数
% 定义6种颜色（MATLAB默认颜色顺序，区分度好）
colors = [0, 0.4470, 0.7410;       % 蓝色
          0.8500, 0.3250, 0.0980;  % 橙红色
          0.9290, 0.6940, 0.1250;  % 黄色
          0.4940, 0.1840, 0.5560;  % 紫色
          0.4660, 0.6740, 0.1880;  % 绿色
          0.6350, 0.0780, 0.1840]; % 深红色

% 定义6种线型（纯线型，不带标记）
linestyles = {'-', '--', ':', '-.', '--', ':'};

% 定义6种标记样式（尺寸适中，避免过重）
% markers = {'o', 's', '^', 'd', 'v', '>'};
marker_sizes = [2, 1, 2, 1, 2, 1];  % 控制标记大小
marker_edge_width = 1;  % 标记边缘线宽
 
    lt = 6;   
    gamm = 13;  %正则项系数
    params.gamm = gamm; 
    max_iter = 240;
    params.max_iter = max_iter ;
    % 预分配结果数组
    tau = [1,1,0,1,0,2];% GPSB1,GPSB2,RPSB,DA,PSM,SM
    theta =  [1,1,1,0.5,1,2];%t_value_1234  = tau.*k.^(theta); 
    k_num = length(tau);
    theta_0 = 0.8*ones(1,k_num);
    params.num = k_num;
    k_rpsb =zeros(lt,k_num); 
    gapx_rpsb = zeros(max_iter,lt,k_num); 
    Fv_rpsb = zeros(max_iter,lt,k_num); 
    gapFF_rpsb = zeros(max_iter,lt,k_num); 
    cpu_time_rpsb = zeros(max_iter,lt,k_num);
    k_time_rpsb = zeros(max_iter,lt,k_num);
    valueF_rpsb = zeros(lt,k_num);  
    cpu_time_sum = zeros(lt,k_num);
    tol = 0.01; 
    params.tol = tol; 
    params.theta_0 = theta_0;%区间为(0,1)
    params.tau = tau;
    params.theta = theta; 
    gammaall = [];

    % 
    % pv = [20,200,400,800, 1000, 2000, 3000, 4000, 6000]; % 对应论文里的特征维度
    % nv = [10,50,200,80, 400, 500, 1200, 2000, 3000];   % 对应论文里的样本数
for t = 7:9
    if t == 4
        gamma = [8e1,8e1,1e1,2e1,2e1,1.5e1]; 
        % gamma = 10*ones(1,6);
        gammaall = [gammaall;gamma];
       elseif t == 5
        gamma = [4e1,6e1,1e1,2e1,1e1,1.5e1]; 
        % gamma = 1e1*ones(1,6);
        gammaall = [gammaall;gamma];
    elseif t == 6
        gamma = [4e1,6e1,1e1,2e1,1e1,1.5e1];  
        % gamma = 1e1*ones(1,6);
        gammaall = [gammaall;gamma];
    elseif t == 7
        gamma = [6e1,8e1,1e1,2e1,1e1,1.5e1];  
        % gamma = 1e1*ones(1,6);
        gammaall = [gammaall;gamma];
    elseif t == 8
        gamma = [6e1,8e1,1e1,2e1,1e1,1.5e1]; 
        % gamma = 1e1*ones(1,6);
        gammaall = [gammaall;gamma];
    elseif t == 9
        gamma = [6e1,8e1,1e1,2e1,1e1,1.2e1]; 
        % gamma = 1e1*ones(1,6);
        gammaall = [gammaall;gamma];
    end
    params.gamma = gamma; 
    i = t-3;
    data = data2_Friedman_Fixed(t);
    A = data.A;
    b = data.b;
    original_signal = data.x;
    m = data.p;
    n = data.n; 

    params.original_signal = original_signal;
    params.m = m; 
    params.A = A;
    params.b = b; 
    params.n = n;            
    rng(1);
    params.x0 = zeros(m,1);
    x_aver = zeros(m,lt,k_num);  
    params.Ab = A'*b;   
    params.nb = 1/(2 * n) * (- norm(b)^2);
    params.Sigmma_h = -(1/n) * (A' * A); 

            % 定义函数句柄
            % $f(x)=\gamma\|x\|_1+\frac{1}{n} {\|Ax-b\|_1}; \\h(x)=-\frac{1}{2n}  \|Ax - b\|^2_2$
            f_subgrad = @(x) HT_f_subgrad(x,params);
            h_prox = @(v,lam) HT_h_prox(v,lam); 
            F = @(x) - 0.5 * (1/n) * norm(A*x-b)^2 + gamm *norm(x,1) + (1/n)*norm(A*x-b, 1);
%             f_prox = @(t,x0,lambda_prox,k,s_k) ...
%                 HT_f_prox(A,x0,lambda_prox,k,s_k,params,t);   
            h_subgrad = @(x) HT_h_subgrad(x, params);
            f_prox_Sigmma = @(v,k,lambda_prox,ki) ...
                HT_f_sigmmaf_prox(v, k,lambda_prox,ki, params);
            for ki = [1,2,4,6]
                params.ki = ki;
                [k_rpsb(i,ki),Fv_rpsb(:,i,ki),gapFF_rpsb(:,i,ki),gapx_rpsb(:,i,ki),x_aver(:,i,ki), cpu_time_rpsb(:,i,ki),valueF_rpsb(i,ki)] = ...
                    HT_PSB(F,f_subgrad,h_subgrad,f_prox_Sigmma,params);
                    cpu_time_sum(i,ki) = sum(cpu_time_rpsb(:,i,ki)); 
            end 
        
  
    % % 创建并设置gap图形 gapx(k, :)  = norm(x-x_true,2);
    % figGap = figure(3*t-11);
    % clf(figGap);
    % set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
    % hold on;legend_label={'GPSB1','GPSB2','DA','SM'};
    % for ki = [1,2,4,6] % 假设有k_num种画法
    %     plot(2:5:k_rpsb(i,ki),gapx_rpsb(2:5:k_rpsb(i,ki),i,ki),...
    %         'Marker',markers{ki},...
    %         'LineStyle', linestyles{ki}, ...
    %         'Color', colors(ki, :), ... 
    %         'LineWidth', 1.5); 
    %     %根据caseki生成图例标签
    %     method_names  = {'GPSB1','GPSB2','DA','SM'};
    % end
    % hold off;    
    % 
    % % xlim([2,k_rpsb(i,ki)-1]);
    % % ylim([0, 5.1]);  
    % legend(method_names, 'Interpreter', 'latex', 'FontSize', 20);
    % xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 18);
    % ylabel('$\|x_k-x_{true}\|$','Interpreter','latex', 'FontSize', 18);  
    % title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold');
    % grid on;
    % set(gca, 'FontSize', 16);   % 增大刻度标签字体，例如 14
    % set(gcf, 'Position', [50, 50, 900, 700]);

    %创建并设置gap图形 gapFF(k,:) = norm(F(history(k, :)')-F(history(k-1, :)')); 
    figGap = figure(3*t-10);
    clf(figGap);
    set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
    hold on;
    for ki = [1,2,4,6] % 假设有k_num种画法 
        plot(1:1:(k_rpsb(i,ki)-1),gapFF_rpsb(1:1:(k_rpsb(i,ki)-1),t-3,ki),... 
            'LineStyle', linestyles{ki}, ...
            'Color', colors(ki, :), ... 
            'LineWidth', 1.5);

        %根据caseki生成图例标签
        method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
        legend_label  = {'GPSB1','GPSB2','DA','SM'};
    end
    hold off;
    %xlim([2,k_rpsb(i,ki)-1]);
    % ylim([0, 1]);  
    legend(legend_label, 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('$k$','Interpreter','latex', 'FontSize', 12);
    ylabel('$\|F(x_k)-F(x_{k-1})\|$','Interpreter','latex', 'FontSize', 12);  
    title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
     
%   %创建并设置gap图形  Fv = F(x_old)
    figGap = figure(3*t-9);
    clf(figGap);
    set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
    hold on;
    for ki = [1,2,4,6] % 假设有k_num种画法
        plot(2:5:(k_rpsb(i,ki)-1),log(Fv_rpsb(2:5:(k_rpsb(i,ki)-1),i,ki)),...  
            'LineStyle', linestyles{ki}, ...
            'Color', colors(ki, :), ... 
            'LineWidth', 2.5);

        %根据caseki生成图例标签
        method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
        legend_label  = {'GPSB1','GPSB2','DA','SM'};
    end
    hold off;
%     xlim([2,k_rpsb(i,ki)-1]);ylim([0, 0.1]);  
    legend(legend_label, 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 18);
    ylabel('$F(x_k)$','Interpreter','latex', 'FontSize', 18);  
    title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 16);   % 增大刻度标签字体 
    ax = gca; 
    ax.YAxis.Exponent = 3;
    set(gcf, 'Position', [850, 150, 900, 700]);


end

%%
%  %创建柱状图 六种维度，四个算法，运行时间，cpu_time_sum 是6*4矩阵
%  figure(1)
% dimensions = [100, 200, 300, 400,500,600]; % 维度
% % datas = cpu_time_sum(:,[1,2,4,6]);
% datas = zeros(6,6);
% for ii=1
%     for j=[1,2,4,6]
%         datas(ii,j) = gapFF_rpsb(k_rpsb(ii,j),ii,j);
%     end
% end 
% datas = datas(:,[1,2,4,6]);
% b = bar(dimensions,datas, 'grouped');
% 
% % 设置颜色（使用MATLAB默认的颜色，可自定义）
% colors = [0, 0.4470, 0.7410;       % 蓝色
%           0.8500, 0.3250, 0.0980;  % 橙红色
%           0.9290, 0.6940, 0.1250;  % 黄色
%           0.4940, 0.1840, 0.5560;  % 紫色
%           0.4660, 0.6740, 0.1880;  % 绿色
%           0.6350, 0.0780, 0.1840]; % 深红色
% 
% for i = 1:length(b)
%     b(i).FaceColor = colors(i,:);
%     b(i).FaceAlpha = 0.8; % 设置透明度
% end
% 
% % 添加标签和标题
% xlabel('Dimension of A', 'FontSize', 12, 'Interpreter','latex');
% ylabel('Error Tolerance $\epsilon_k$', 'FontSize', 12, 'Interpreter','latex');
% % title('$\beta_k = 2\sqrt{k}$','Interpreter','latex', 'FontSize', 14);
% 
% % 设置纵坐标范围（根据您的数据调整）
% % ylim([0,max(max(datas))/3+1]);
% 
% % 添加图例
% % legend('$(\theta_0,tau,\theta)=(0.1,0.8,0.1)$','$(\theta_0,tau,\theta)=(0.3,0.8,0.3)$', '$(\theta_0,tau,\theta)=(0.5,0.8,0.5)$', '$(\theta_0,tau,\theta)=(0.7,0.8,0.7)$', ...
% %        'Location', 'northeast', 'FontSize', 12, 'Interpreter','latex');
% legend( 'GPSB1', 'GPSB2','DA','SM', ...
%        'Location', 'northeast', 'FontSize', 12, 'Interpreter','latex');
% 
% % 设置坐标轴属性
% dimensions_labels = { '800*80', '1000*400', '2000*500', '3000*1200', '4000*2000', '6000*3000'}; % 问题维度
% set(gca, 'XTick', 100:100:100*length(dimensions), 'XTickLabel', dimensions_labels);
% set(gca, 'FontSize', 12, 'LineWidth', 1);
% grid on;
% set(gcf, 'Position', [100, 100, 800, 600]);
% 


%%
%  %创建柱状图 六种维度，四个算法,迭代次数，k_rpsb是6*4矩阵
% 
%  figure(2)
%  datas2 = k_rpsb(:,[1,2,4,6]);
%  dimensions = [100, 200, 300, 400,500,600]; % 维度
% b = bar(dimensions, datas2, 'grouped');
% 
% % 设置颜色
% colors = [0, 0.4470, 0.7410;       % 蓝色
%           0.8500, 0.3250, 0.0980;  % 橙红色
%           0.9290, 0.6940, 0.1250;  % 黄色
%           0.4940, 0.1840, 0.5560;  % 紫色
%           0.4660, 0.6740, 0.1880;  % 绿色
%           0.6350, 0.0780, 0.1840]; % 深红色
% 
% for i = 1:length(b)
%     b(i).FaceColor = colors(i,:);
%     b(i).FaceAlpha = 0.8; % 设置透明度
% end
% 
% % 添加标签和标题
% xlabel('Dimension of A', 'FontSize', 12, 'Interpreter','latex');
% ylabel('Iteration k', 'FontSize', 12, 'Interpreter','latex');
% % title('$\beta_k = 2\sqrt{k}$','Interpreter','latex', 'FontSize', 14);
% 
% % 设置纵坐标范围
% % ylim([0,max(max(k_rpsb))/3+10]);
% 
% % 添加图例 
% legend('DA','RPSB', 'GPSB1', 'GPSB2', ...
%        'Location', 'northeast', 'FontSize', 12, 'Interpreter','latex');
% 
% % 设置坐标轴属性
% dimensions_labels = { '800*80', '1000*400', '2000*500', '3000*1200', '4000*2000', '6000*3000'}; % 问题维度
% set(gca, 'XTick', 100:100:100*length(dimensions), 'XTickLabel', dimensions_labels);
% set(gca, 'FontSize', 12, 'LineWidth', 1);
% grid on;
% set(gcf, 'Position', [100, 100, 800, 600]);


 %%
%     %单独画线性图    存图 存数据
%     mv = [20,200,400,800,1000,2000,3000,4000,6000];
%     nv = [10,50,120,240,400,700,1000,1200,2800];
% for t = 4:9
%     i = t-3;
% 
%         % 创建并设置gap图形
%         k_time_rpsb = zeros(200,6,6);
%         for i0 = 1:k_rpsb
%             k_time_rpsb(i0) = sum(1:i0);
%         end
%         m = mv(t);
%         n = nv(t);
%         figGap = figure(i+6);
%         clf(figGap);
%         set(figGap,'Name',['d=',num2str(mv(t)),' n=',num2str(nv(t))]);
%         hold on;
%         for ki = [1,2,4,6] % 假设有k_num种画法
%             plot(5:5:k_rpsb(i,ki),Fv_rpsb(5:5:k_rpsb(i,ki),i,ki),...
%                 linestyles{ki},'Color',colors(ki,:),'MarkerSize',2,'LineWidth',2); 
%             % 根据caseki生成图例标签
%             method_names = {'DA','RPSB', 'GPSB1',  'GPSB2'};
%         legend_label  = {'GPSB1','GPSB2','DA','SM'};
%         end
%         hold off;
%     %     xlim([2,20*i]);
%     %     ylim([0, 0.8]);  
%         legend(legend_label, 'Interpreter', 'latex', 'FontSize', 18);
%         xlabel('$k$','Interpreter','latex', 'FontSize', 18);
%         ylabel('$F(x_k)$','Interpreter','latex', 'FontSize', 18);  
%         title(['$(m,n) = $','(',num2str(mv(t)),',',num2str(nv(t)),')'],'Interpreter','latex', 'FontSize', 18, 'FontWeight', 'bold');
%         grid on;
% 
% %     %创建并设置gap图形 gapFF(k,:) = norm(F(history(k, :)')-F(history(k-1, :)'));
% %     figGap = figure(3*t-10);
% %     clf(figGap);
% %     set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
% %     hold on;
% %     for ki = [1,2,4,6] % 假设有k_num种画法
% %         plot(1:1:(k_rpsb(i,ki)-1),gapFF_rpsb(1:1:(k_rpsb(i,ki)-1),t-3,ki),... 
% %             'LineStyle', linestyles{ki}, ...
% %             'Color', colors(ki, :), ... 
% %             'LineWidth', 1.5);
% %     
% %         %根据caseki生成图例标签
% %         method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
% %         legend_label  = {'GPSB1','GPSB2','DA','SM'};
% %     end
% %     hold off;
% %     %xlim([2,k_rpsb(i,ki)-1]);
% %     ylim([0, 1]);  
% %     legend(legend_label, 'Interpreter', 'latex', 'FontSize', 12);
% %     xlabel('$k$','Interpreter','latex', 'FontSize', 12);
% %     ylabel('$\|F(x_k)-F(x_{k-1})\|$','Interpreter','latex', 'FontSize', 12);  
% %     title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 12, 'FontWeight', 'bold');
% %     grid on;
% % 
% % %   %创建并设置gap图形  Fv = F(x_old)
% %     figGap = figure(2*i-1);
% %     clf(figGap);
% %     set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
% %     hold on;
% %     for ki = [1,2,4,6] % 假设有k_num种画法
% %         plot(2:10:(k_rpsb(i,ki)-1),Fv_rpsb(2:10:(k_rpsb(i,ki)-1),i,ki),...
% %             'LineStyle', linestyles{ki}, ...
% %             'Color', colors(ki, :), ... 
% %             'LineWidth', 2.5);
% % 
% %         %根据caseki生成图例标签
% %         method_names = {'GPSB1','GPSB2','RPSB','DA','PSM','SM'};
% %         legend_label  = {'GPSB1','GPSB2','DA','SM'};
% %     end
% %     hold off;
% % %     xlim([2,k_rpsb(i,ki)-1]);ylim([0, 0.1]);  
% %     legend(legend_label, 'Interpreter', 'latex', 'FontSize', 20);
% %     xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 18);
% %     ylabel('$F(x_k)$','Interpreter','latex', 'FontSize', 18);  
% %     title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold');
% %     grid on;
%     set(gca, 'FontSize', 16);   % 增大刻度标签字体 
%     ax = gca; 
% 
% %     ax.XAxis.Exponent = [];
%     ax.YAxis.Exponent = 3;
%     set(gcf, 'Position', [850, 150, 800, 600]);
% % 
% %     % 创建并设置gap图形 gapx(k, :)  = norm(x-x_true,2);
% %     figGap = figure(2*i);
% %     clf(figGap);
% %     set(figGap,'Name',['m=',num2str(m),' n=',num2str(n)]);
% %     hold on;legend_label={'GPSB1','GPSB2','DA','SM'};
% %     for ki = [1,2,4,6] % 假设有k_num种画法
% %         plot(2:10:k_rpsb(i,ki),gapx_rpsb(2:10:k_rpsb(i,ki),i,ki),...
% %             'Marker',markers{ki},...
% %             'LineStyle', linestyles{ki}, ...
% %             'Color', colors(ki, :), ... 
% %             'LineWidth', 1.5); 
% %         %根据caseki生成图例标签
% %         method_names  = {'GPSB1','GPSB2','DA','SM'};
% %     end
% %     hold off;    
% % 
% %     % xlim([2,k_rpsb(i,ki)-1]);
% %     % ylim([0, 5.1]);  
% %     legend(method_names, 'Interpreter', 'latex', 'FontSize', 20);
% %     xlabel('Iteration $k$','Interpreter','latex', 'FontSize', 18);
% %     ylabel('$\|x_k-x_{true}\|$','Interpreter','latex', 'FontSize', 18);  
% %     title(['$(m,n) = $','(',num2str(m),',',num2str(n),')'],'Interpreter','latex', 'FontSize', 20, 'FontWeight', 'bold');
% %     grid on;
% %     set(gca, 'FontSize', 16);   % 增大刻度标签字体，例如 14
% %     set(gcf, 'Position', [50, 50, 800, 600]);
% % 
% % 
% 
% end

%%
    % 存图 
         save_folder = 'D:\BaiduSyncdisk\3 weakly convex PSB\MATLABcodes\MATLABcodeHT_review1\20260629\';  %实际路径
        % 确保文件夹存在，如果不存在则创建
         if ~exist(save_folder,'dir')
             mkdir(save_folder);
         end 

    figures = findall(groot,'Type','figure');
    for idx = 1:length(figures) 
        % 切换到当前图形
        figure(figures(idx));    
        % 获取标题（尝试从坐标轴获取）
        fileName = [num2str(idx),'_4Psb6DimSameIterFv']; 
        % 处理非法字符
        fileName = regexprep(fileName,'[\\/:*?"<>|]','_'); 

        % 完整文件路径
        fullPath = fullfile(save_folder,fileName);   
    %     保存为EPS格式（矢量图）
        print(figures(idx),[fullPath '.eps'],'-depsc','-r600');   
        % 保存为FIG格式
        savefig(figures(idx),[fullPath '.fig']);
    end
    % 存数据到单个文件mat 
%             定义目标文件夹路径（示例：保存到 D 盘的 'MATLAB_Data' 文件夹）
         save_folder = 'D:\BaiduSyncdisk\3 weakly convex PSB\MATLABcodes\MATLABcodeHT_review1\20260629\';  % 请替换为你的实际路径
        % 确保文件夹存在，如果不存在则创建
         if ~exist(save_folder,'dir')
             mkdir(save_folder);
         end 
%      拼接完整文件路径并保存
         save_file = fullfile(save_folder,'ALL_4Psb6DimSameIter.mat');
         save(save_file,"colors","cpu_time_sum","cpu_time_rpsb","F","h_subgrad","h_prox",...
             "valueF_rpsb","f_prox_Sigmma","f_subgrad","gamm","gammaall","Fv_rpsb","gapx_rpsb","valueF_rpsb",...
             "gapFF_rpsb","k_num","k_rpsb","legend_label","linestyles","method_names","tol",...
             "markers","theta_0","theta","tau");  % 保存所有工作区变量

