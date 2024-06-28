%For x^(-s) in [1e-8,1]

function [B,XT,Leb,Sg] = REIM_power2(M)
    a = 1e-8; b = 1;
    c = 1e-5; d = 1e-2;
    X =unique([a:(c-a)/3000:c,c:(d-c)/3000:d,d:(b-d)/3000:b]');
    
    N = length(X);
    g = @(b) 1./(X+b);
    E = 10.^(linspace(-10,1,1000))';
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