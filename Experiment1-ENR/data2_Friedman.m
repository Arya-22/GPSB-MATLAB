function data = data2_Friedman(opts)
% 'Regularization Paths for Generalized Linear Models via Coordinate Descent'
rng(40); 
n = 720;
p = 2560;
sigma = 0.5;
rho = 0.5;
Hchoice = 1;  %  1:'Identity', 2:'Orthogonal'
 
Sigmachoice = 2; %  1:'CS', 2:'AR(1)'
 

% 读取参数
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'n');         n = opts.n;                       end
if isfield(opts, 'p');         p = opts.p;                       end
if isfield(opts, 'sigma');     sigma = opts.sigma;               end
if isfield(opts, 'rho');       rho = opts.rho;                   end
if isfield(opts, 'Hchoice');   Hchoice= opts.Hchoice;            end
if isfield(opts, 'Sigmachoice');   Sigmachoice= opts.Sigmachoice;            end


% 协方差矩阵的结构
if Sigmachoice == 1
    Sigma = rho + (1-rho)*eye(p); %compound symmetric structures (CS)
elseif Sigmachoice == 2
    Sigma = zeros(p,p);
    for i=1:p
        for j=1:p
            Sigma(i,j)=rho^(abs(i-j));
        end
    end
end

% 生成数据
A = mvnrnd(zeros(p,1),Sigma,n); % or A = mvnrnd(zeros(p,1),eye(p),n)*chol(Sigma);
A = zscore(A,0,1); % 对矩阵A进行标准化，使其每一列的均值为0，标准差为1。
Hx=zeros(p,1);
for i=1:p
    Hx(i) = (-1)^(i)*exp(-(2*i-1)/20);
end
if Hchoice == 1
    H = eye(p);
    x = Hx;
elseif Hchoice == 2
    H=randn(p,p);
    [H,~,~]=svd(H);
    x = H'*Hx;
end

err = sigma*randn(n,1);
b=A*x+err;

if sigma~=0
    SNR = var(A*x)/var(err);
else
    SNR = 0;
end
data.m=p;
data.n=n;
data.A=A;
data.x=x;
data.H=H;
data.b=b;
data.SNR = SNR;
end