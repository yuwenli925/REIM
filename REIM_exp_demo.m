% For exp(-tau x) and varphi(-tau x) in [lambda_min,lambda_max]

% phi = @(z) (exp(z) - 1)./z;
phi = @(z) (exp(z) - 1)./z;

lambda_min = 1; lambda_max = 1e6;
tau = linspace(1e-3,1,1e2);

M = 30;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max,'exp');

Xtest = linspace(lambda_min,lambda_max,5e5)';
gtest = 1./(Xtest+Bm');
errphi = zeros(length(tau),1); 
errexp = errphi;
for i = 1:length(tau)
    phiz = phi(-tau(i)*Xtest);
    phizi = phi(-tau(i)*Xm);
    exz = exp(-tau(i)*Xtest);
    exzi = exp(-tau(i)*Xm);
    errphi(i) = norm(phiz - gtest*(Gm\phizi), 'inf');
    errexp(i) = norm(exz - gtest*(Gm\exzi), 'inf');
end
semilogy(tau,errexp,'o','MarkerSize',5,'Color',[0.00 0.45 0.74],'LineWidth',1)
fprintf('The maximum error for exp(-x) is %e\n',max(errexp));
hold on
semilogy(tau,errphi,'r*','MarkerSize',5,'LineWidth',1)
fprintf('The maximum error for phi(-x) is %e\n',max(errphi));
xlabel('$\tau_i$','interpreter','latex','fontsize',16)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
legend('$\exp(-\tau_i x)$','$\varphi(-\tau_i x)$','interpreter','latex','fontsize',12)
