lambda_min = 1e-6; lambda_max = 1;
M = 30;
[B,XT,L,Sg] = eim(M,lambda_min,lambda_max);
figure(1)
semilogy(1:M,sort(Sg),'o','Color',[0.00 0.45 0.74]);
xlabel('$i$','interpreter','latex','fontsize',16)
ylabel('$b_i$','interpreter','latex','fontsize',16)
figure(2)
semilogy(1:M,sort(XT),'r*');
xlabel('$i$','interpreter','latex','fontsize',16)
ylabel('$x_i$','interpreter','latex','fontsize',16)
figure(3)
plot(1:M,L,'-square','Color',[0.93 0.69 0.13],'LineWidth',1);
xlabel('$n$','interpreter','latex','fontsize',16)
ylabel('$L_n$','interpreter','latex','fontsize',16)

function [B,XT,L,Sg] = eim(M,a,b)
    X = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/4000:b]');
    N = length(X);
    g = @(b) 1./(X+b);
    E = 10.^(linspace(-7,1,1000))'; 
    Q = zeros(N,M); %q
    Sg = zeros(M,1); %u1,u2...uM
    Wg = zeros(N,M); %g1,g2.
    T = 0; %tM
    B = 0; %BijM
    Sg(1) = E(1); %- poles
    Wg(:,1) = g(Sg(1)); %g1
    [~,T(1)] = max(abs(Wg(:,1))); 
    Q(:,1) = Wg(:,1)/Wg(T(1),1); %q1
    B = Q(T(1),1);
    L(1) = max(abs(Q(:,1)/B)); %Lebesgue constants
    for m = 2:M
        l = zeros(length(E),1);
        for i = 1:length(E) %argmax res
            gtrain = g(E(i,:));
            fi = Q(T,1:m-1)\gtrain(T);
            l(i) = max(abs(gtrain - Q(:,1:m-1)*fi));
        end
        [~,index] = max(l);
        Sg(m) = E(index);%um
        Wg(:,m) = g(Sg(m,:));
        K = Wg(T,m);
        Sig = B\K;
        rM = Wg(:,m) - Q(:,1:(m-1))*Sig;
        [~,T(m)] = max(abs(rM));
        Q(:,m) = rM/rM(T(m));
        B = Q(T,1:m);
        L(m) = norm(Q(:,1:m)/B,"inf");
    end
    XT = X(T); %interpolation points
end