function Rate_single_user = Single_Rate_Comp(Gamma,BETAA,Pu,N,Phi,C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[M,K] = size(Gamma);

PhiPhi = zeros(K,K);
PC_2 = zeros(K,K);

for ii=1:K
     for k=1:K
         PhiPhi(ii,k) = round( Phi(:,k)'*Phi(:,ii) );
     end
end

for m_iter = 1:M
    A = zeros(M,K);
    A(m_iter,:) = ones(1,K);
    Gamma_new = Gamma;
    for ii=1:K
         for k=1:K
             PC_2(ii,k) = sum( (C(ii)*A(m_iter,k)*Gamma_new(m_iter,k)/BETAA(m_iter,k))*BETAA(m_iter,ii)  )*PhiPhi(ii,k);
         end
    end
    PC22=N^2*(abs(PC_2)).^2;


    for k=1:K
        deno1=0;
        for m=1:M
            deno1=deno1 + C(k)*A(m,k)*Gamma_new(m,k)*sum(BETAA(m,:));
        end
        IN(k) = (N*sum(Gamma_new(:,k).*A(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
        

        if C(k)==1
            Rate_single_user(m_iter,k) = log2( 1+Pu*N^2*( sum(C(k)*A(:,k).*Gamma_new(:,k)) )^2/IN(k) );
        end
    end 

end
end

