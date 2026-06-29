%定义次梯度
function grad = HT_f_subgrad(x,params)
%    $f(x)=\gamma\|x\|_1+\frac{1}{n} {\|Ax-b\|_1}; \\h(x)=-\frac{1}{2n}  \|Ax - b\|^2_2$
    %计算f(ω)的次梯度 
    A = params.A;
    b = params.b;
    [n,~] = size(A); 
    % 梯度: -\frac{1}{n}A(Ax-b)
        y = A*x - b;
        u =  - A' * y ./n;
    grad =  u;
end
