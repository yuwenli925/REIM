% For ratioanl approximation of x^(-s) in [1e-6,1] and [1e-8,1]
% xm records the interpolation points, bm records the negative pole

function [xm,bm,G] = REIM(M,a,b)
    if a == 1e-6 && b == 1
        xset = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/4000:b]');
        bset = 10.^(linspace(-7,1,1000))';
    elseif a == 1e-8 && b == 1
        c = 1e-5; d = 1e-2;
        xset =unique([a:(c-a)/3000:c,c:(d-c)/3000:d,d:(b-d)/3000:b]');
        bset = 10.^(linspace(-10,1,1000))';
    else 
        error("The sets of x and b should be redefined.");
    end
    % xset = linspace(a,b,1e4)';
    N = length(xset);
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
