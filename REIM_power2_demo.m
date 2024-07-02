%% For x^(-s) in [lambda_min,lambda_max]

s = 0.5;
lambda_min = 1e-8; lambda_max = 1;

M = 31;
[Xm,Bm,Gm] = REIM(M,lambda_min,lambda_max);

Xtest = linspace(lambda_min,lambda_max,5e5)';
gtest = 1./(Xtest+Bm');

err = norm(Xtest.^(-s) - gtest*(Gm\Xm.^(-s)), 'inf');

fprintf('The error is %e\n',max(err));

function [xm,bm,G] = REIM(M,a,b)
    c = 1e-5; d = 1e-2;
    xset =unique([a:(c-a)/3000:c,c:(d-c)/3000:d,d:(b-d)/3000:b]');
    % xset = linspace(a,b,1e4)';
    N = length(xset);
    bset = 10.^(linspace(-10,1,1000))';
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
%     [bm, id] = sort(bm);
%     xm = xm(id);
%     G = G(id,id);
end