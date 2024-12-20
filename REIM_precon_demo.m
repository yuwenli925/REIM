%For (x^(-0.5)+Kx^0.5)^-1 in [1e-6,1]
K = sort(rand(100,1)*(1-1e-6)+1e-6);
M = 30; s = 0.5;
[Xm,Bm,Gm] = REIM(M,1e-6,1,"precon");
X = linspace(1e-6,1,5e5)';
errREIM =zeros(100,1); 
gx = 1./(repmat(Xm,1,length(Bm))+Bm');
gtest = 1./(repmat(X,1,length(Bm))+Bm');
for i = 1:length(K)
  fx = 1./(Xm.^(-s)+K(i)*Xm.^s);
  ftest = 1./(X.^(-s)+K(i)*X.^s);
  res = gx\fx;
  errREIM(i) = max(abs(ftest - gtest*res));
end
semilogy(K,errREIM,'o')
xlabel('$K_i$','interpreter','latex','fontsize',16)
ylabel('$L_{\infty}$ error','interpreter','latex','fontsize',16)

