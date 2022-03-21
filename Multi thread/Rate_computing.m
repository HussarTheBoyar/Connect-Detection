function [rate_user,A,status,num,sort_order,reject] = Rate_computing(Gamma,BETAA,Phi,Pu,RReq_co,Req,N)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);
rate = zeros(M,K);
rate_user = zeros(1,K);
PC_2 = zeros(K,K);
IN   = zeros(1,K);
A    = zeros(M,K);
status = zeros(1,K);
sort_order     = zeros(M,K);    %sort the addresses of APs with gamma from large to small
reject = zeros(1,K);            %the unconsidered order in the CVX condition string
frac   = zeros(1,K);            %gamma_check / gamma_req ratio, through sorting this value we will get reject

for ii=1:K
    for k=1:K
        PC_2(ii,k) = sum( (Gamma(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phi(:,k)'*Phi(:,ii);
    end
end
PC22=N^2*(abs(PC_2)).^2;

if K == 1
    PC22 = 0;
end

for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + Gamma(m,k)*sum(BETAA(m,:));
    end
    IN(k) = (sum(Gamma(:,k)) + Pu*deno1*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
    for m = 1:M
        rate(m,k) = Pu*( Gamma(m,k) )^2/IN(k);
    end
end

%computing Gamma require
RReq = RReq_co*Req;
Gamma_req = zeros(1,K);
for k = 1:K
    SINR_req = 2^RReq-1;
    Gamma_req(k) = sqrt(SINR_req*IN(k)/Pu);
end    

for k =1:K
    [tmp_Gamma, sort_order(:,k)] = sort(Gamma(:,k));  
    Gamma_check = 0;
    
    for m = 1:M
        Gamma_check = Gamma_check+tmp_Gamma(m);
        if Gamma_check>=Gamma_req(k)
            A( sort_order(m,k),k ) = 1;
            break
        else
            A( sort_order(m,k),k ) = 1;
        end
    end
    
    if Gamma_check<Gamma_req(k)
        status(k) = 1;
    end
    
    frac(k) = Gamma_check/Gamma_req(k);
    
end

for k = 1:K
    Gamma_new = Gamma(:,k).*A(:,k);
    rate_user(k) = N^2*Pu*sum( Gamma_new )^2/IN(k);
end

%Number of APs connect to user
num = sum(A);

%Get element to reject
[~,reject] = sort(frac);

end

