%For x^(-s) in [1e-6,1]
s = 0.25; %s = 0.5;0.75;0.95;
M = 40;
[~,XT,Leb,Sg] = REIM(M);
fx = XT.^(-s);
gx = 1./(repmat(XT,1,length(Sg))+Sg');
X = linspace(1e-6,1,5e5)';
ftest = X.^(-s);
gtest = 1./(repmat(X,1,length(Sg))+Sg');
err = []';
for i = 1:M    
 res = gx(1:i,1:i)\fx(1:i);
 err(i) = max(abs(ftest - gtest(:,1:i)*res));
end
semilogy(1:M,err)
hold on
xlabel('$n$','interpreter','latex','fontsize',18)
ylabel('$L^{\infty}$ error','interpreter','latex','fontsize',18)

function [B,XT,Leb,Sg] = REIM(M)
    a = 1e-6; b = 1;
    X = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/4000:b]');
    N = length(X);
    g = @(b) 1./(X+b);
    E = 10.^(linspace(-7,1,1000))';
    Q = zeros(N,M); %q
    Sg = zeros(M,1); %-1*poles
    Wg = zeros(N,M); %g1,g2...gM
    T = 0; %tM
    B = 0; %BijM
    Sg(1) = E(1); %u1
    Wg(:,1) = g(Sg(1)); %g1
    [~,T(1)] = max(abs(Wg(:,1))); 
    Q(:,1) = Wg(:,1)/Wg(T(1),1); %q1
    B = Q(T(1),1);
    Leb(1) = max(abs(Q(:,1)/B));%Lebesgue constant
    for m = 2:M
        L = zeros(length(E),1);
        for i = 1:length(E) %argmax res
            fi = Q(:,1:m-1)\g(E(i,:));
            L(i) = max(abs(g(E(i,:)) - Q(:,1:m-1)*fi));
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
    XT = X(T);%interpolation points
end