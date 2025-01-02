% For x^(-s) in [lambda_min,lambda_max] = [1e-6,1]

s = [0.25,0.5,0.75,0.95];
lambda_min = 1e-6; lambda_max = 1;
M = 40;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max,'power');
X = unique([lambda_min:(0.001-lambda_min)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(lambda_max-0.01)/4000:lambda_max]');
Xtest = linspace(lambda_min,lambda_max,5e5)';
err = zeros(M,1);
for j = 1:length(s)
    fx = X.^(-s(j));
    ftest = Xtest.^(-s(j));
    %%rEIM
    for i = 1:M
        gtest = 1./(Xtest+Bm(1:i)');
        err(i) = norm(ftest - gtest*(Gm(1:i,1:i)\Xm(1:i).^(-s(j))), 'inf');
    end
    fprintf('The rEIM error for s = %.2f is %e\n',s(j),err(end));
    figure(j)
    semilogy(1:M,err,'-o','LineWidth',1)
    title(sprintf("$s = %.2f$", s(j)),'interpreter','latex','fontsize',16)
    xlabel('$n$','interpreter','latex','fontsize',16)
    ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
    hold on
    %%OGA
    [~, ~, errOGA] = OGA(s(j),30);
    semilogy(1:30,errOGA,'-*','LineWidth',1)
    hold on 
    %%AAA
    f = @(x) x.^(-s(j));
    errAAA2 = []';
    for i = 1:M
        [~,pole,~,~,~,~,~,errAAA] = aaa(f,X,mmax=i);
        if (i >=2 && errAAA(end) > 10*errAAA(end-1))||length(errAAA) < i %quit AAA
            iter = i-1;
            break
        end
        gx = 1./(repmat(X,1,length(pole))-pole');
        gtest = 1./(repmat(Xtest,1,length(pole))-pole');
        errAAA2(i) = max(abs(ftest-gtest*(gx\fx)));
        iter = min(i,errAAA);
    end
    semilogy(1:iter,errAAA(1:iter),'-square','LineWidth',1)
    hold on
    semilogy(1:iter,errAAA2(1:iter),'-+','LineWidth',1)
    legend('rEIM','OGA','AAA','AAA*','interpreter','latex','fontsize',12)
end

