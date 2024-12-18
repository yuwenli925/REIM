% For 1/(x^s+d) in [1e-6,1]

s = [0.5, 1];

lambda_min = 1e-6; lambda_max = 1;
d = sort(rand(100,1)*(1e3-1)+1); 
Lambda = 1e6;

M = 30;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max);

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
fprintf('The error for s = 0.5 is %e\n',max(err));
for i = 1:length(d)
    phiz = 1./(Xtest.^s(2)+d(i)/Lambda^s(2));
    phizi = 1./(Xm.^s(2)+d(i)/Lambda^s(2));
    err(i) = norm(phiz - gtest*(Gm\phizi), 'inf');
end
figure(2)
semilogy(d,err,'r*','MarkerSize',5)
fprintf('The error for s = 1 is %e\n',max(err));
xlabel('$d_i$','interpreter','latex','fontsize',16)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',16)
function [xm,bm,G] = REIM(M,a,b)
    c = b/1e4; d = b/1e2;
    xset = unique([a:(c-a)/2000:c,c:(d-c)/2000:d,d:(b-d)/3000:b]');
    % xset = linspace(a,b,1e4)';
    N = length(xset);
    % bset = 10.^(linspace(-3,3,1000))'; 
    bset = 10.^(linspace(-6,2,1000))';
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