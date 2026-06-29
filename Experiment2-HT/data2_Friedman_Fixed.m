function data = data2_Friedman_Fixed(case_id, opts)
    % 设置随机种子以保证结果可复现
    rng(42); 
    
    pv = [20,200,400,800, 1000, 2000, 3000, 4000, 6000]; % 对应你论文里的特征维度
    nv = [10,50,200,80, 400, 500, 1200, 2000, 3000];   % 对应你论文里的样本数
    cases = [pv', nv'];
    
    p = cases(case_id, 1);
    n = cases(case_id, 2);
    fprintf('生成第 %d 种情况数据: p= %d, n = %d\n', case_id, p, n);
    
    % --- 默认参数 ---
    sigma = 0.5;   % 噪声缩放系数
    rho = 0.5;     % 特征相关性系数
    Sigmachoice = 2; % 1:'CS', 2:'AR(1)'
    noise_type = 't3'; % 'gaussian', 't3' (Student-t df=3), 'cauchy'
    sparsity_ratio = 0.05; % 稀疏度，例如 5%
    
    % 读取自定义参数
    if ~exist('opts', 'var'); opts = []; end
    if isfield(opts, 'n'); n = opts.n; end
    if isfield(opts, 'p'); p = opts.p; end
    if isfield(opts, 'sigma'); sigma = opts.sigma; end
    if isfield(opts, 'rho'); rho = opts.rho; end
    if isfield(opts, 'Sigmachoice'); Sigmachoice = opts.Sigmachoice; end
    if isfield(opts, 'noise_type'); noise_type = opts.noise_type; end
    
    % --- 1. 生成设计矩阵 A (保留你的优秀做法) ---
    if Sigmachoice == 1
        Sigma = rho + (1-rho)*eye(p); % CS
    elseif Sigmachoice == 2
        Sigma = zeros(p, p); % AR(1)
        for i=1:p
            for j=1:p
                Sigma(i,j) = rho^(abs(i-j));
            end
        end
    end
    A = mvnrnd(zeros(p,1), Sigma, n); 
    A = zscore(A, 0, 1); % 标准化
    
    % --- 2. 生成严格稀疏的 Ground Truth 信号 x (解决审稿人痛点) ---
    s = max(5, round(p * sparsity_ratio)); % 非零元素个数
    x_true = zeros(p, 1);
    nonzero_indices = randperm(p, s);
    x_true(nonzero_indices) = randn(s, 1); % 非零元素服从正态分布
    % 为了使得真实信号更强，也可以 x_true(nonzero_indices) = sign(randn(s,1)) + randn(s,1);
    
    % --- 3. 生成真正的重尾噪声 err (解决审稿人痛点) ---
    if strcmp(noise_type, 'gaussian')
        err = sigma * randn(n, 1); % 仅供对比试验使用
    elseif strcmp(noise_type, 't3')
        err = sigma * trnd(3, n, 1); % 自由度为3的t分布，典型的重尾噪声
    elseif strcmp(noise_type, 'cauchy')
        err = sigma * trnd(1, n, 1); % 柯西分布，极端重尾
    else
        error('未知的噪声类型');
    end
    
    % 生成观测值 b
    b = A * x_true + err;
    
    % 计算 SNR (对于重尾分布，方差可能无穷大，这里算伪SNR作参考)
    if sigma ~= 0
        SNR = var(A * x_true) / var(err(abs(err) < prctile(abs(err),95))); % 剔除极值算近似SNR
    else
        SNR = 0;
    end
    
    % 封装返回数据
    data.n = n;
    data.p = p;
    data.A = A;         % 设计矩阵 (对应你公式里的 \Phi)
    data.x = x_true; % 严格稀疏真实信号
    data.b = b;         % 观测值
    data.SNR = SNR; 
end