% For 1/(x^s+d) in [1e-6,1]

s = [0.5, 1];

lambda_min = 1e-6; lambda_max = 1;
d = sort(rand(100,1)*(1e3-1)+1); 
Lambda = 1e6;

M = 30;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max,'time');

Xtest = linspace(lambda_min,lambda_max,5e5)';
gtest = 1./(Xtest+Bm');
err = zeros(length(d),1); 
for i = 1:length(d)
    phiz = 1./(Xtest.^s(1)+d(i)/Lambda^s(1));
    phizi = 1./(Xm.^s(1)+d(i)/Lambda^s(1));
    err(i) = norm(phiz - gtest*(Gm\phizi), 'inf');
end
figure(1)
semilogy(d,err,'o','MarkerSize',5,'Color',[0.00 0.45 0.74])
xlabel('$d_i$','interpreter','latex','fontsize',16)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
fprintf('The maximum error for s = 0.5 is %e\n',max(err));
for i = 1:length(d)
    phiz = 1./(Xtest.^s(2)+d(i)/Lambda^s(2));
    phizi = 1./(Xm.^s(2)+d(i)/Lambda^s(2));
    err(i) = norm(phiz - gtest*(Gm\phizi), 'inf');
end
figure(2)
semilogy(d,err,'r*','MarkerSize',5)
fprintf('The maximum error for s = 1 is %e\n',max(err));
xlabel('$d_i$','interpreter','latex','fontsize',16)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
