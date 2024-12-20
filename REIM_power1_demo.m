% For x^(-s) in [lambda_min,lambda_max] = [1e-6,1]

s = [0.25,0.5,0.75,0.95];
lambda_min = 1e-6; lambda_max = 1;
M = 40;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max,'power');
Xtest = linspace(lambda_min,lambda_max,5e5)';
err = zeros(M,1);
for j = 1:4
    for i = 1:M
        gtest = 1./(Xtest+Bm(1:i)');
        err(i) = norm(Xtest.^(-s(j)) - gtest*(Gm(1:i,1:i)\Xm(1:i).^(-s(j))), 'inf');
    end
    fprintf('The error for s = %.2f is %e\n',s(j),err(end));
    figure(j)
    semilogy(1:M,err,'-o','Color',[0.00 0.45 0.74],'LineWidth',1)
    title(sprintf("$s = %.2f$", s(j)),'interpreter','latex','fontsize',16)
    xlabel('$n$','interpreter','latex','fontsize',16)
    ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
end

