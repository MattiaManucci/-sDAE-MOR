%% Code based on the work:
%[1] M. Manucci and B. Unger, Balancing-based model reduction for switched descriptor systems
% ArXiv e-print 2404.10511, 2024.
%% --------------------------------------------------------------
function [S,rho] = Solve_LS_GLE(A,S,T,B,nf,tol,flag)
tol_tSV=0.1*tol; %Truncation tollerance
opts.maxit=1e3; opts.tol=1e-4; %SVD convergence options
%% Solve for Reachability Gramian
if flag==1
    n=size(B,1); 
    %% Computing Constant beta for the rescaling
    JfunInv=@(v,flag) JfunInv2(A{1},S{1},T{1},nf{1},v,flag); %Function for inv(J1)*v and inv(J1)'*v
    sigma_min=svds(JfunInv,[n,n],1,'largest',opts); %This computes the largest SV of A^-1
    sigma_min=1/sigma_min;
    normD=0;

    for i=1:numel(A)
        Dfun3=@(v,flag) Dfun2(A{i},A{1},S{i},S{1},T{i},T{1},nf{i},nf{1},v,flag);
        normD=normD+svds(Dfun3,[n,n],1,'largest',opts)^2;
    end
    rho=(2*normD/(2*sigma_min));
    %% --------------------------------------------------------------
    tol2=tol*sigma_min/3; % Stopping Criteria for Krylov Subspace (KS)
    nit=200; % Max number of iterations for KS
    %If norm of known term is larger than one then rescale with the
    %constant rho.
    normB=normest(B);
    if normB>1
        B=sqrt(1/rho)*B;
    end
    [Vs,Y] = solve_KS(A{1},S{1},T{1},B,nf{1},tol2,nit);
    [U,Sig,~]=svd(Y,'econ');
    [s,~]=find(diag(Sig)<tol_tSV);
    if isempty(s)
        s=size(Y,2)+1;
    end
    Z_l=Vs*U(:,1:(s-1))*sqrt(Sig(1:(s-1),1:(s-1)));
    maxit=20; % Maximum number of iterations for GLEs solver
    %% Main Loop for GLE solver
    for j=1:maxit
        B_n=[];
        for i=2:numel(A)
            DP_l=Dfun(A{i}, A{1}, S{i}, S{1}, T{i}, T{1}, nf{i}, nf{1}, Z_l);
            B_n=[B_n, DP_l];
        end
        % Rescale the known term
        B_n=[B,(1/sqrt(rho))*B_n];

        [Vs,Y] = solve_KS(A{1},S{1},T{1},B_n,nf{1},tol2,nit);
        [U,Sig,~]=svd(Y,'econ'); 
        [s,~]=find((diag(Sig))<tol_tSV); % Extract the relevant information from Y

        if isempty(s)
            s=size(Y,2)+1;
        end
        
        Z_l_new=Vs*U(:,1:(s-1))*sqrt(Sig(1:(s-1),1:(s-1)));
       
        if size(Z_l_new,2)>size(Z_l,2)
            diff=size(Z_l_new,2)-size(Z_l,2);
            Z_l=[Z_l,zeros(size(Z_l,1),diff)];
        else
            Z_l=Z_l(:,1:size(Z_l_new,2));
        end
        GLE_ET=((norm(Z_l,'fro')+norm(Z_l_new,'fro'))*norm(abs(Z_l_new)-abs(Z_l),'fro'));
        display(GLE_ET)
        % Exit tolerance for the method
        if GLE_ET<0.5*tol
            
            fprintf('Solution of generalized Lypunov for Reachability Gramian converged \n')
            S=Z_l_new;
            break
        end
        Z_l=Z_l_new;
    end
    if j>=maxit
        
        fprintf('Solution of Generalized Lypunov for Reachability Gramian did not converge, Error Estimate is: %d \n',GLE_ET)
        S=Z_l_new;
    end
end
%% Solve for Observability Gramian
if flag==2

    n=size(B,2); % size of J1
    %% Computing constant beta for the rescaling
    JfunInv=@(v,flag) JfunInv3(A{1},S{1},T{1},nf{1},v,flag); %function for inv(J1')*v and inv(J1)*v
    sigma_min=svds(JfunInv,[n,n],1,'largest',opts); %This actually compute the largest SV of J1'^-1
    sigma_min=1/sigma_min;
    normD=0;
    for i=1:numel(A)
        Dfun4=@(v,flag) Dfun5(A{i},A{1},S{i},S{1},T{i},T{1},nf{i},nf{1},v,flag);
        normD=normD+svds(Dfun4,[n,n],1,'largest',opts)^2;
    end
    rho=(2*normD/(2*sigma_min));
    %% --------------------------------------------------------------
    tol2=tol*sigma_min/3; % Stopping Criteria for Krylov Subspace (KS)
    nit=200; % Max number of iterations for KS
    normB=normest(B');
    %If norm of known term is larger than one then rescale with the
    %constant rho.
    if normB>1
        B=sqrt(1/rho)*B;
    end

    [Vs,Y] = solve_KS_t(A{1},S{1},T{1},B',nf{1},tol2,nit);
    [U,Sig,~]=svd(Y,'econ'); 
    [s,~]=find((diag(Sig))<tol_tSV); % Extract the relevant information from Y
 
    if isempty(s)
        s=size(Y,2)+1;
    end
    S_l=Vs*U(:,1:(s-1))*sqrt(Sig(1:(s-1),1:(s-1))); 
    maxit=20; 
    for j=1:maxit
        B_n=[];
        for i=2:numel(A)
            DP_l=Dfun_t(A{i}, A{1}, S{i}, S{1}, T{i}, T{1}, nf{i}, nf{1}, S_l);
            B_n=[B_n, DP_l];
        end
        B_n=[B',sqrt(1/rho)*B_n];

        [Vs,Y] = solve_KS_t(A{1},S{1},T{1},B_n,nf{1},tol2,nit);
        [U,Sig,~]=svd(Y,'econ'); [s,~]=find((diag(Sig))<tol_tSV);
     
        if isempty(s)
            s=size(Y,2)+1;
        end
        if s(1)==1
            break
        end
        S_l_new=Vs*U(:,1:(s-1))*sqrt(Sig(1:(s-1),1:(s-1)));

        if size(S_l_new,2)>size(S_l,2)
            diff=size(S_l_new,2)-size(S_l,2);
            S_l=[S_l,zeros(size(S_l,1),diff)];
        else
            S_l=S_l(:,1:size(S_l_new,2));
        end

        % Error Estimate
        GLE_ET=((norm(S_l,'fro')+norm(S_l_new,'fro'))*norm(abs(S_l_new)-abs(S_l),'fro'));
        display(GLE_ET)
        %% Exit criteria
        if GLE_ET<0.5*tol
            
            fprintf('Solution of generalized Lypunov for Observability Gramian converged \n')
            S=S_l_new;
            break
        end
        S_l=S_l_new;
    end
    if j>=maxit
        
        fprintf('Solution of Generalized Lypunov for Observability Gramian did not converge, error estimate is: %d \n',GLE_ET)
        S=S_l_new;
    end
end
end

%% Functions for Matrix-Vecotr multiplication: A*v and A'*v

function [Av] = Jfun(A, S, T, nf,v)
    nc=size(v,2);
    v=[v;zeros(size(A,1)-nf,nc)];
    Tv=T*v; SATv=S\(A*(Tv)); 
    Av=SATv(1:nf,:);
end

function [Avt] = Jfun_t(A, S, T, nf, v)
    nc=size(v,2);
    v=[v;zeros(size(A,1)-nf,nc)];
    Stv=S'\v; TtAtStv=T'*A'*Stv;
    Avt=TtAtStv(1:nf,:);
end

function [Av] = JfunInv(A, S, T, nf,v)
    nc=size(v,2);
    v=[v;zeros(size(A,1)-nf,nc)];
    Sv=S*v; ATv=T\(A\(Sv)); 
    Av=ATv(1:nf,:);
end

function [Avt] = JfunInv_t(A, S, T, nf, v)
    nc=size(v,2);
    v=[v;zeros(size(A,1)-nf,nc)];
    Ttv=T'\v; Avt=S'*(A'\Ttv);
    Avt=Avt(1:nf,:);
end


function [v]=JfunInv2(A, S, T, nf, v, flag)
     switch flag
         case 'notransp'
             v=JfunInv(A, S, T, nf, v);
         case 'transp'
             v=  JfunInv_t(A, S, T, nf, v);
     end
end


function [v]=JfunInv3(A, S, T, nf, v, flag)
     switch flag
         case 'notransp'
             v=JfunInv_t(A, S, T, nf, v);
         case 'transp'
             v=  JfunInv(A, S, T, nf, v);
     end
end
%% Functions for Matrix-Vecotr multiplication: D*v and D'*v

function [Dv] = Dfun(A, A1, S, S1, T, T1, nf, nf1, v)

    nc=size(v,2);

    J1v=Jfun(A1, S1, T1, nf1, v);
    
    v=[v;zeros(size(A,1)-nf1,nc)];

    Pv=T\(T1*v); Pv=Pv(1:nf,:);
    Pv=[Pv;zeros(size(A,1)-nf,nc)];

    Jiv=S\(A*T*Pv); Jiv=Jiv(1:nf,:);
    Jiv=[Jiv;zeros(size(A,1)-nf,nc)];

    Jiv=T1\(T*Jiv); Jiv=Jiv(1:nf1,:);

    Dv=Jiv-J1v;

end

function [Dvt] = Dfun_t(A, A1, S, S1, T, T1, nf, nf1, v)
    
    nc=size(v,2);
    J1vt=Jfun_t(A1, S1, T1, nf1, v);

    v=[v;zeros(size(A,1)-nf1,nc)];

    Pvt=T'*(T1'\v); Pvt=Pvt(1:nf,:);
    Pvt=[Pvt;zeros(size(A,1)-nf,nc)];

    Jivt=T'*A'*(S'\Pvt); Jivt=Jivt(1:nf,:);
    Jivt=[Jivt;zeros(size(A,1)-nf,nc)];

    Jivt=T1'*(T'\Jivt); Jivt=Jivt(1:nf1,:);
  
    Dvt=Jivt-J1vt;

end

function [v]=Dfun2(A,A1,S,S1,T,T1,nf,nf1,v,flag)
     switch flag
         case 'notransp'
             v=  Dfun(A, A1, S, S1, T, T1, nf, nf1, v);
         case 'transp'
             v=Dfun_t(A, A1, S, S1, T, T1, nf, nf1, v);
     end
end

function [v]=Dfun5(A,A1,S,S1,T,T1,nf,nf1,v,flag)
     switch flag
         case 'notransp'
             v=    Dfun_t(A, A1, S, S1, T, T1, nf, nf1, v);
         case 'transp'
             v=      Dfun(A, A1, S, S1, T, T1, nf, nf1, v);
     end
end