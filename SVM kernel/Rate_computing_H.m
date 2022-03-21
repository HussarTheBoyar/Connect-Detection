function [rate_user,status,num,RReq] = Rate_computing_H(Gamma,BETAA,Phi,Pu,RReq_co,Req,N)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
rate_user = zeros(1,K);
PC_2 = zeros(K,K);
IN   = zeros(1,K);
A    = ones(M,K);
status = zeros(1,K);
for ii=1:K
    for k=1:K
        PC_2(ii,k) = sum( (Gamma(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phi(:,k)'*Phi(:,ii);
    end
end
PC22=(abs(PC_2)).^2*N^2;


for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + Gamma(m,k)*sum(BETAA(m,:));
    end
    IN(k) = (sum(Gamma(:,k))*N + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
end

RReq = RReq_co*log2(1+Req);

for k = 1:K
    Gamma_new = Gamma(:,k);
    rate_user(k) = log2(1+ Pu*N^2*sum( Gamma_new )^2/IN(k));
    if rate_user(k) < RReq
        status(k) = 1;
    end
end

%Number of APs connect to user
num = sum(A);


end

