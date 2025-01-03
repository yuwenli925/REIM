%% Solve space-fractional parabolic PDE by rational approximation and adaptive time step size

addpath './FEM'
[Xm,Bm,~] = REIM(30,1e-6,1,'time');

s = 0.5; %%change s = 0.5 or 1 to generate Fig.6
Lambda = 1e6; 
gx = 1./(repmat(Xm,1,length(Bm))+Bm');
tend = 1; tol = 1e-4;
tau = 0.001; 
T = 0; tau0 = tau;
np = length(Bm);
[node, elem] = squaremesh([0, 1, 0, 1],0.125);
for iter=1:5
    [node, elem] = uniformrefine(node, elem);
end
NV = size(node,1);
aux = myauxstructure(elem);
[D, area] = gradbasis(node,elem);
[Stiff, Mass] = P1mat2d(D,area,elem);
bv = unique(aux.bdEdge);
fv = setdiff(1:NV,bv)';
S = Stiff(fv,fv);
M = Mass(fv,fv);

LL = cell(np,1); UU = LL; PP = LL; QQ = LL;
e=ones(np,1);
for i=1:np
    [LL{i}, UU{i}, PP{i}, QQ{i}] = lu(S/Lambda + Bm(i)*M);
end

err = zeros(length(T),1);
Err2 =err;
errest = err;
Tdel = []; taudel=[]; 

Uarray = cell(length(T),1);
uexact = Uarray; fn = Uarray;
uexact{1} = u(0,node(fv,:));
Uarray{1} = uexact{1};


while T(end) <= tend
    j = length(T)+1;
    uexact{j} = u(T(end)+tau0,node(fv,:));
    rhs = P1rhs2d(node, elem, area, @(x0) f(T(end)+tau0,x0,s));
    fn{j} = rhs(fv);

    fx = 1./(Xm.^s+1/(tau0*Lambda^s));
    res = gx\fx;
    Uj = Uarray{j-1};
    F = ((M*Uj)/tau0 + fn{j})/Lambda^s;
    %BDF2
    a = tau(end);
    k0=(a+2*tau0)/(tau0*(a+tau0)); 
    k1=-(a+tau0)/(a*tau0);
    k2=tau0/(a*(a+tau0));
    if j>2
        Ui = Uarray{j-2};
        G = (-k1*(M*Uj)-k2*(M*Ui) + fn{j})/Lambda^s;
    else 
        G = F; 
    end
    hx = 1./(Xm.^s+k0/Lambda^s);
    res2 = gx\hx; 

    U1 = zeros(length(fv),1); U2 = U1;
    for i = 1:np
        solve = QQ{i}*( UU{i}\( LL{i}\(PP{i}*[F G]) ) );
        U1 = U1 + res(i)*solve(:,1);
        U2 = U2 + res2(i)*solve(:,2);
    end
    if j ==2, U2=U1; end
    %errest(j-1) = max(abs(U1-U2));
    errest(j-1) = sqrt((U1-U2)'*M*(U1-U2));
    if errest(j-1) <= tol
        tau(j-1) = tau0; 
        T(j) = T(j-1) + tau0;
        Uarray{j} = U1; 
        %err(j) = max(abs(U1-uexact{j}));
        err(j) = sqrt((U1-uexact{j})'*M*(U1-uexact{j}));
        tau0 = 0.8*tau0*(tol/max(errest(j-1),1e-6))^0.5;     
       if T(j)+tau0 > tend               
               if T(j-1) >= tend
                   break 
               end
               tau0 = tend - T(j);               
        end  
    else 
        Tdel = [Tdel;T(j-1)+tau0];
        taudel = [taudel;tau0];
        tau0 = 0.8*tau0*(tol/errest(j-1))^0.5;
    end
    if tau0 <= 1e-4
           break
    end
    %T(end)
end
  
figure(1)
plot(T,err,'-','Color',[0.00 0.45 0.74],'LineWidth',1)
hold on
plot(T(1:4:end),err(1:4:end),'o','Color',[0.00 0.45 0.74],'LineWidth',1)
xlim([0.001,1])
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$L^2$ error','interpreter','latex','fontsize',18)

figure(2)
plot(T(1:4:end),tau(1:4:end),'square-','Color',[0.00 0.45 0.74],'LineWidth',1)
hold on
xlim([0,1])
scatter(Tdel,taudel,'x','LineWidth',1)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$\tau_n$','interpreter','latex','fontsize',18)
legend('Accepted Step Sizes','rejected Step Sizes','interpreter','latex','fontsize',12)

%--------------------------------------------------------------------------
function z = u(t,p)
z = exp(-t/20)*cos(2*pi*t)*sin(pi*p(:,1)).*sin(pi*p(:,2));
end

function z = f(t,p,s)
z = (-1/20+(2*pi^2)^s)*u(t,p)-2*pi*exp(-t/20)*sin(2*pi*t)*sin(pi*p(:,1)).*sin(pi*p(:,2));
end

