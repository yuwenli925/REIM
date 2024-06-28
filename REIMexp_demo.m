%For exp(-taux) and phi(-taux) in [0,1e6]
tau = sort(rand(100,1)*(1-1e-3)+1e-3);
M = 30;
[~,XT,~,Sg] = REIM(M);
gx = 1./(repmat(XT,1,length(Sg))+Sg');
X = linspace(0,1e6,5e5)';
gtest = 1./(repmat(X,1,length(Sg))+Sg');
err = []';
for i = 1:length(tau)
    fx = exp(-tau(i)*XT); %exp(-taux)
    ftest = exp(-tau(i)*X);
    % fx = (1-exp(-tau(i)*XT))./(tau(i)*XT);fx(isnan(fx)) = 1; %phi(-taux)
    % ftest = (1-exp(-tau(i)*X))./(tau(i)*X);ftest(isnan(ftest)) = 1;
    res = gx\fx;
    err(i) = max(abs(ftest - gtest*(res)));
end
semilogy(tau,err,'o')
hold on

function [B,XT,Leb,Sg] = REIM(M)
    a = 0; b = 1e6; c=1e2; d = 1e4;
    X = unique([a:(c-a)/2000:c,c:(d-c)/2000:d,d:(b-d)/3000:b]');
    N = length(X);
    g = @(b) 1./(X+b);
    E = 10.^(linspace(-3,4,1000))'; 
    Q = zeros(N,M); %q
    Sg = zeros(M,1); %u1,u2...uM
    Wg = zeros(N,M); %g1,g2...gM
    T = 0; %tM
    B = 0; %BijM
    Sg(1) = E(1); %u1
    Wg(:,1) = g(Sg(1)); %g1
    [~,T(1)] = max(abs(Wg(:,1))); 
    Q(:,1) = Wg(:,1)/Wg(T(1),1); %q1
    B = Q(T(1),1);
    Leb(1) = max(abs(Q(:,1)/B));
    for m = 2:M
        L = zeros(length(E),1);
        for i = 1:length(E) %argmax res
            % fi = lsqminnorm(Q(:,1:m-1),g(E(i,:)));
            fi = Q(:,1:m-1)\g(E(i,:));
            L(i) = max(abs(g(E(i,:)) - Q(:,1:m-1)*fi));
            %L(i) = Lp(g(E(i,:)) - Q(:,1:m-1)*fi,p);
        end
        [~,index] = max(L);
        Sg(m) = E(index);%um
        Wg(:,m) = g(Sg(m,:));
        K = Wg(T,m);
        Sig = B\K;
        rM = Wg(:,m) - Q(:,1:(m-1))*Sig;
        [~,T(m)] = max(abs(rM));
        Q(:,m) = rM/rM(T(m));
        B = Q(T,1:m);
        Leb(m) = norm(Q(:,1:m)/B,"inf");
    end
    XT = X(T);
end