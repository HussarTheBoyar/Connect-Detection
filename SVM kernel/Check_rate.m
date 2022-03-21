function [rate_user] = Check_rate(Gamma,BETAA,Phi,Pu,N,A,C)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
rate_user = zeros(1,K);
PC_2 = zeros(K,K);


for ii=1:K
    for k=1:K
        PC_2(ii,k) = sum( (C(ii)*A(:,k).*Gamma(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phi(:,k)'*Phi(:,ii);
    end
end
PC22=N^2*(abs(PC_2)).^2;

if K == 1
    PC22 = 0;
end

for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + C(k)*A(m,k)*Gamma(m,k)*sum(BETAA(m,:));
    end
    IN(k) = (N*sum(Gamma(:,k).*A(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k)); 
end

%computing Gamma require

for k = 1:K
    Gamma_new = Gamma(:,k).*A(:,k);
    if C(k) == 1
        rate_user(k) = log2(1+ Pu*N^2*sum( Gamma_new )^2/IN(k));
    end
end

end

