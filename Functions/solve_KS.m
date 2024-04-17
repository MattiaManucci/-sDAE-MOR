%% Code based on the work:
%[1] M. Manucci and B. Unger, Balancing-based model reduction for switched descriptor systems
% ArXiv e-print 2404.10511, 2024.
%% --------------------------------------------------------------
function [V,Y] = solve_KS(A,S,T,B,nf,tol,k) % All this function has a cost that is O(n)
% Iteration 0
n=size(A,1);
Vk=B; V=[];
nGS=3; tolGS=1e-6;  I=speye(n);
for ii=1:size(Vk,2)
    for j=1:nGS
        if (ii>1)&&(isempty(V)~=1)
            Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
            Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
        end
        Vk(abs(Vk(:,ii))<1e-12,ii)=0;
    end
    if norm(Vk(:,ii))>tolGS
        Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
        V=[V,Vk(:,ii)];
    end
end
I_n1=I(1:nf,:); I_n1(:,(nf+1):end)=0;
V0=V; V=V0; n=size(V,2); AV=I_n1*(S\(A*T*I_n1'*V)); VAV=V'*AV; BV=V'*B;
for i=1:k

    Vk=I_n1*(S\(A*T*I_n1'*V0));
    for ii=1:size(Vk,2)

        for j=1:nGS
            Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
            Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
            Vk(abs(Vk(:,ii))<1e-12,ii)=0;
        end

        if norm(Vk(:,ii))>tolGS
            Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
            V=[V,Vk(:,ii)];
        end
    end
    n_new=size(V,2);
    if i>1
        Vj=V(:,n+1:n_new);
        Res=2*Vj'*(AV*Y);
        KS_RES=norm(Res,"fro");
        if KS_RES<tol
            V=V(:,1:n);
            break
        end
    end

    V0=V(:,n+1:n_new);  AVnew=I_n1*(S\(A*T*I_n1'*V0));

    VAV=[       VAV,        V(:,1:n)'*AVnew ;...
        V(:,n+1:n_new)'*AV ,V(:,n+1:n_new)'*AVnew];

    AV=[AV,AVnew];
    BV=[BV;V0'*B];
    Y=lyap(VAV,BV*BV');
    n=n_new;
end
end