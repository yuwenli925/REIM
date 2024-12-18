%For (x^(-0.5)+Kx^0.5)^-1 in [1e-6,1]
K = sort(rand(100,1)*(1-1e-6)+1e-6);
M = 30; s = 0.5;
[Xm,Bm,Gm] = REIM(M,1e-6,1);
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
semilogy(K,errREIM,'o','LineWidth',1)
xlabel('$K_i$','interpreter','latex','fontsize',16)
ylabel('$L_{\infty}$ error','interpreter','latex','fontsize',16)

function [xm,bm,G] = REIM(M,a,b)
    xset = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/4000:b]');
    % xset = linspace(a,b,1e4)';
    N = length(xset);
    bset = 10.^(linspace(-6,2,1000))';
    % bset = 10.^(linspace(0,3,1000))';
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
end