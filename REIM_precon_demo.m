%For (x^(-0.5)+Kx^0.5)^-1 in [1e-6,1]
K = sort(rand(100,1)*(1-1e-6)+1e-6);
M = 30;
[~,XT,~,Sg] = REIM(M);
X = linspace(1e-6,1,5e5)';
errREIM =zeros(100,1); 
gx = 1./(repmat(XT,1,length(Sg))+Sg');
gtest = 1./(repmat(X,1,length(Sg))+Sg');
for i = 1:length(K)
  fx = 1./(XT.^(-0.5)+K(i)*XT.^0.5);
  ftest = 1./(X.^(-0.5)+K(i)*X.^0.5);
  res = gx\fx;
  errREIM(i) = max(abs(ftest - gtest*res));
end
semilogy(d,errREIM,'o')

function [B,XT,Leb,Sg] = REIM(M)
    a = 1e-6; b = 1;
    X = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/4000:b]');
    N = length(X);
    g = @(b) 1./(X+b);
    E = 10.^(linspace(-6,2,1000))'; 
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