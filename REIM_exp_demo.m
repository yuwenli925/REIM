% For exp(-tau x) and varphi(-tau x) in [lambda_min,lambda_max]

% phi = @(z) (exp(z) - 1)./z;
phi = @(z) (exp(z) - 1)./z;

lambda_min = 1; lambda_max = 1e6;
tau = linspace(1e-3,1,1e2);

M = 30;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max);

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
fprintf('The error is %e\n',max(err));
hold on
semilogy(tau,errphi,'r*','MarkerSize',5,'LineWidth',1)
xlabel('$\tau_i$','interpreter','latex','fontsize',16)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
legend('$\exp(-\tau_i x)$','$\varphi(-\tau_i x)$','interpreter','latex','fontsize',12)
function [xm,bm,G] = REIM(M,a,b)
    c = b/1e4; d = b/1e2;
    xset = unique([a:(c-a)/2000:c,c:(d-c)/2000:d,d:(b-d)/3000:b]');
    % xset = linspace(a,b,1e4)';
    N = length(xset);
    % bset = 10.^(linspace(-3,3,1000))'; 
    bset = 10.^(linspace(0,3,1000))';
    bm = zeros(M,1); %u1,u2...uM
    xm = zeros(M,1);
    gm = zeros(N,M); %g1,g2...gM
    bm(1) = bset(1); %u1
    gset = 1./(xset+bset');
    gm(:,1) = gset(:,1);
    [~,id] = max(abs(gm(:,1))); 
    xm(1) = xset(id);
    G = 1/(xm(1)+bm(1));
    for m = 2:M
        L = zeros(length(bset),1);
        for i = 1:length(bset) % argmax res
            L(i) = norm( gset(:,i) - gm(:,1:m-1)*(G\(1./(xm(1:m-1)+bset(i)))) , 'inf');
        end
        [~,id] = max(L);
        bm(m) = bset(id); % um
        gm(:,m) = gset(:,id);
        rm = gm(:,m) - gm(:,1:m-1)*(G\(1./(xm(1:m-1)+bm(m))));
        [~,id] = max(abs(rm));
        xm(m) = xset(id);
        G = 1./(xm(1:m)+bm(1:m)');
    end
    [bm, id] = sort(bm);
    xm = xm(id);
    G = G(id,id);
end
