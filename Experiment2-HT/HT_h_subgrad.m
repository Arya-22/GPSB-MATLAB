function subgrad = HT_h_subgrad(x,params)
    % $f(x)=\gamma\|x\|_1+\frac{1}{n} {\|Ax-b\|_1}; \\h(x)=-\frac{1}{2n}  \|Ax - b\|^2_2$
    %计算f(ω)的次梯度 
    
    A = params.A;
    b = params.b;
    gamm = params.gamm;
    omega = x;
    axb = A*x - b;
    [n,~] = size(A); 
    %     % L1范数的次梯度
    l1_subgrad1 = zeros(length(omega), 1);
    for i = 1:length(omega)
        if omega(i) > 0
            l1_subgrad1(i) = 1;
        elseif omega(i) < 0
            l1_subgrad1(i) = -1;
        else
            % 在零点，次梯度是[-1, 1]中的任意值
            l1_subgrad1(i) = 2*rand(1,1)-1; % 可以取0或[-1,1]中的其他值
        end
    end
    
    %     % L1范数的次梯度
    l1_subgrad2 = zeros(length(axb), 1);
    for i = 1:length(axb)
        if axb(i) > 0
            l1_subgrad2(i) = 1;
        elseif axb(i) < 0
            l1_subgrad2(i) = -1;
        else
            % 在零点，次梯度是[-1, 1]中的任意值
            l1_subgrad2(i) = 2*rand(1,1)-1; % 可以取0或[-1,1]中的其他值
        end
    end
    subgrad = gamm*l1_subgrad1 + 1/n * A'*l1_subgrad2;
   
end



