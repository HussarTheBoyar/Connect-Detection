function [rate_user,A,status,num,RReq] = Rate_computing_TWC(Gamma,BETAA,Phi,Pu,RReq_co,Req,N)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
rate_user = zeros(1,K);
PC_2 = zeros(K,K);
IN   = zeros(1,K);
A    = zeros(M,K);
rp = zeros(1,K);

H = sqrt(Gamma);
H_check = mean(H);

for k = 1:K
    for m = 1:M
        if H(m,k)>=H_check(k)
            A(m,k) = 1;
        end
    end
end

Gamma = Gamma.*A;

for ii=1:K
    for k=1:K
        PC_2(ii,k) = sum( (A(:,ii).*Gamma(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phi(:,k)'*Phi(:,ii);
    end
end
PC22=N^2*(abs(PC_2)).^2;

if K == 1
    PC22 = 0;
end

for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + A(m,k)*Gamma(m,k)*sum(BETAA(m,:));
    end
    IN(k) = N*(sum(Gamma(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k)); 
end

%computing Gamma require
RReq = RReq_co*log2(1+Req);



for k = 1:K
    Gamma_new = Gamma(:,k).*A(:,k);
    rate_user(k) = log2(1+ Pu*N^2*sum( Gamma_new )^2/IN(k));
end

for k=1:K
    if rate_user(k)<RReq
        rp(k) = 1;
        rate_user(k) = 0;
    end
end


%Number of APs connect to user
num = sum(A);
status = sum(rp)/K;
end

