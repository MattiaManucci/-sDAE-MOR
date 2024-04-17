%% Compute reduced solution and errors
rv=input('Size of the reduced problems: \n'); %You can also inset a vector here
nr=numel(rv);
ns=input('Insert number of switches \n');
ns=ns+1;
seed = 123; rng(seed); % For the first switching path
%seed = 1234;rng(seed); % For the Second switching path
% Choose random order of switches
rr = 1:ns; rr=mod(rr,nm+1); rr(rr==0)=nm;
random_indices = randperm(length(rr));
rr = rr(random_indices); 
% Time scale and time step
scale=0.1; dt=scale*1e-1;
%% Inputs
% First input signal
u=@(t) kron(ones(size(B1,2),1),sin(t.^2+t));
du=@(t) kron(ones(size(B1,2),1), cos(t.^2+t).*(2*t+1));
ddu=@(t) kron(ones(size(B1,2),1),-sin(t.^2+t).*(2*t+1).^2+cos(t.^2+t)*2);
% Second input signal
TT=8;
% u=@(t) kron(ones(size(B1,2),1),sin(2*pi*(exp(t/TT))));
% du=@(t) kron(ones(size(B1,2),1), cos(2*pi*(exp(t/TT))).*(2*pi*exp(t/TT)/TT));
% ddu=@(t) kron(ones(size(B1,2),1), -sin(2*pi*exp(t/TT)).*((2*pi*exp(t/TT))/TT).^2+ cos(2*pi*(exp(t/TT))).*(2*pi*exp(t/TT)/(TT^2)));

%% Third Input Signal (use this signal, scale=1 and  10 switches to reproduce Figure 2 (left) and Figure 1 (right))
% scale=1;
% u=@(t) kron(ones(size(B1,2),1),sin(t));
% du=@(t) kron(ones(size(B1,2),1), cos(t));
% ddu=@(t) kron(ones(size(B1,2),1),-sin(t));
%% ----------------------------------------------------------------------------------------------------
if flag==2
    ddu=@(t) [];
end
%% Compute the solution of the FOM
Dimp=cell(nm,1); dtt=cell(ns,1); t=zeros(ns,1); tspan=cell(ns,1); y=cell(ns,1);
for i=1:nm
    Dimp{i}=-C{i}*[Bi{i}];
end
tic
for i=1:ns
    if i==1
        Q=sparse(n,n);  Q(1:n1{rr(i)},1:n1{rr(i)})=speye(n1{rr(i)});
        x_0=zeros(n,1)-T{rr(i)}*Q*(T{rr(i)}\(Bi{rr(i)}*([u(0);du(0);ddu(0)])));
        t(i)=scale*randi([1 10],1,1);
        dtt{i}=0:dt:t(i);
    else
        t(i)=t(i-1)+scale*randi([1 10],1,1);
        dtt{i}=t(i-1):dt:t(i);
       
        Q=sparse(n,n);  Q(1:n1{rr(i)},1:n1{rr(i)})=speye(n1{rr(i)});
        x_0=(T{rr(i)}*Q*(T{rr(i)}\(x(end,:)'+...
            Bi{rr(i-1)}*([u(tspan{i-1}(end));du(tspan{i-1}(end));ddu(tspan{i-1}(end))])...
            )));
    end
    odefun=@(t,y) DS{rr(i)}*(S{rr(i)}\(A{rr(i)}*y))+Bd{rr(i)}*u(t);

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [tspan{i},x]=ode45(odefun,dtt{i},x_0,options);
    y{i}=Cd{rr(i)}*(x')+Dimp{rr(i)}*(([u(tspan{i}')',du(tspan{i}')',ddu(tspan{i}')'])');
end
time_FOM=toc;

%% Generate the base removing co-linearity for the input-dependent jump
IMP_Jump=[];
for i=1:nm
    for j=1:nm
        if i~=j
            IMP_Jump=[IMP_Jump,T{(i)}*Q*(T{(i)}\Bi{i})];
        end
    end
end
[Ui,Si,Vi]=svds(IMP_Jump,size(IMP_Jump,2)); ind=find(diag(Si)<eps*1e1);
Z_ll=Z_l;
for i=1:3
    Z_ll=Z_ll-Ui(:,i)*(Ui(:,i)'*Z_l);
end

H=Z_l'*S_l; [U,Sss,Vs]=svd(H);
%% Solve ROM
yred=cell(ns,nr);
for jj=1:nr
    r=rv(jj);
    % Projection Matrices VV and W
    VV=Z_l*U(:,1:r)*diag((Ss(1:r)).^(-0.5)); 
    W= S_l*(Vs(:,1:r))*diag((Ss(1:r)).^(-0.5)); % Check right sign for this 
    r=size(VV,2); 
    % Reduced Operatos
    Ared=cell(nm,1); Bred=cell(nm,1); Cred=cell(nm,1); Pired=cell(nm,1);
    for i=1:nm
        n1{i}=size(V{i},2);
        Q=sparse(n,n); Q(1:n1{i},1:n1{i})=speye(n1{i});
        Ared{i}=(W'*VV)\(W'*DS{i}*(S{i}\(A{i}*VV)));
        Bred{i}=(W'*VV)\(W'*Bd{i});
        Cred{i}=Cd{i}*VV;
        Pired{i}=(W'*VV)\(W'*T{i}*Q*(T{i}\VV));
        for j=1:nm
            PiUred{i,j}=(W'*VV)\(W'*T{i}*Q*(T{i}\Bi{j}));
        end
    end
    % Reduced solution size r
    tic
    for i=1:ns

        if i==1
            x_0=zeros(r,1)-PiUred{rr(i),rr(i)}*([u(0);du(0);ddu(0)]);
        else
            x_0=Pired{rr(i)}*(x(end,:)')+PiUred{rr(i),rr(i-1)}*([u(tspan{i-1}(end));du(tspan{i-1}(end));ddu(tspan{i-1}(end))]); 
        end
        odefun=@(t,y) Ared{rr(i)}*y+Bred{rr(i)}*u(t);
        options = odeset('RelTol',1e-10,'AbsTol',1e-10);
        [~,x]=ode45(odefun,dtt{i},x_0,options);
        yred{i,jj}=Cred{rr(i)}*(x')+Dimp{rr(i)}*([u(tspan{i}')',du(tspan{i}')',ddu(tspan{i}')'])';
    end
    toc
end
%% Plot full and prescribed reduced solution
rP=input('Choose the reduced solution for the plot: \n');
Out_to_Plot=input('Insert number of outputs to visualize: \n'); %it can be a single output, a part of the outputs or the all outputs
figure
for i=1:ns

    if i==1
        % Because of the legend here is always plotted the first output
        plot(tspan{i},(y{i}(2,:)'),'--b','LineWidth',LW)
        hold on
        plot(tspan{i},(yred{i,rP}(2,:)'),'-r','LineWidth',LW)

    end
    plot(tspan{i},(y{i}(Out_to_Plot,:)'),'-b','LineWidth',LW)
    hold on

    plot(tspan{i},(yred{i,rP}(Out_to_Plot,:)'),'--r','LineWidth',LW)

end

xlabel('t','Interpreter','Latex')
ylabel('$\mathbf{y}(t)$','Interpreter','Latex')
lgd=legend(['Full Model $N=' num2str(n) '$'],['ROM DAE-SLS $n=' num2str(rv(rP)) '$']);
set(lgd, 'Interpreter','Latex','Location','best');
title(['Number of modes ' num2str(nm) ', Number of switches ' num2str(ns-1)], 'Interpreter', 'Latex')

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');



%% Computing L_2 errors using trapezoidal quadrature rule
final_error_L2_r=zeros(nr,size(Dimp{1},1));
for jj=1:nr
    normU=0; final_error_L2=zeros(1,so); normy=zeros(1,so);
    for i=1:so
        for j=1:ns
            final_error_L2(i)=final_error_L2(i)+(sum(norm(y{j}(:,2:(end-1))-yred{j,jj}(:,2:(end-1))).^2))...
                +0.5*(norm(y{j}(:,1)-yred{j,jj}(:,1)).^2)...
                +0.5*(norm(y{j}(:,end)-yred{j,jj}(:,end)).^2);

            normy(i)=normy(i)+(sum(norm(y{j}(:,2:(end-1))).^2))...
                +0.5*(norm(y{j}(:,1)).^2)...
                +0.5*(norm(y{j}(:,end)).^2);
            if i==1
                normU=normU+sum(u(dtt{j}(2:(end-1))').^2)+0.5*u(dtt{j}(1)').^2+...
                    0.5*u(dtt{j}(end)').^2;
            end
        end
        final_error_L2(i)=sqrt(dt)*sqrt(final_error_L2(i));
        normy(i)= sqrt(dt)*sqrt(normy(i));
    end
    normU=sqrt(dt)*sqrt(normU);
    final_error_L2_r(jj,:)=final_error_L2./normU';
end
%% Generating the plots
figure
semilogy(rv,decay(rv),'-ob','LineWidth',LW)
hold on
semilogy(rv,decay_Sv(rv),'--*r','LineWidth',LW)
semilogy(rv,final_error_L2_r(:,1),'-+k','LineWidth',LW)

xlabel('$r$','Interpreter','Latex')
lgd=legend('$2\sum_{r+1}^{n_{\mathcal{P},\mathcal{Q}}}\tilde{\sigma}_i(\tilde{\mathbf{\mathcal{H}}})$','$2\sum_{r+1}^{n_{\mathcal{P},\mathcal{Q}}}\sigma_i(\tilde{\mathbf{\mathcal{H}}})$','$\frac{\|\mathbf{y}-\tilde{\mathbf{y}}\|_{L_2}}{\|\mathbf{u}\|_{L_2}}$');
set(lgd, 'Interpreter','Latex','Location','best');
set(gca,'Fontname',FN,'Fontsize',1.2*FS);
set(gcf, 'Color', 'w');

