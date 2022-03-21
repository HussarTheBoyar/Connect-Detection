function [PC22,deno1,IN,rate,DS] = PC_comp_new(Gamma,BETAA,Pu,Phi,N,A,C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(Gamma);
PC2 = zeros(K,K);
IN = zeros(1,K);
rate = zeros(1,K);
% A = zeros(M,K);
% C = ones(1,K);
% 
% A(1:2,:) = ones(2,K);

PhiPhi = zeros(K,K);
for k=1:K
    for ii=1:K
        PhiPhi(ii,k) =  Phi(:,k)'*Phi(:,ii);
    end
end

PhiPhi = round(PhiPhi);


for ii=1:K
    for k=1:K
        PC2(ii,k) = sum( ( C(ii).*A(:,k).*Gamma(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*PhiPhi(ii,k);
    end
end

PC22=N^2*(abs(PC2)).^2;

deno1 = zeros(K,1);
for k=1:K
    for m=1:M
        deno1(k)=deno1(k) + C(k)*A(m,k)*Gamma(m,k)*sum(BETAA(m,:));
    end
    IN(k) = (N*sum(Gamma(:,k).*A(:,k)) + Pu*deno1(k)*N + Pu*sum(PC22(:,k)) - Pu*PC22(k,k));
    DS(k) = N^2*Pu*C(k)*sum(A(:,k).*Gamma(:,k))^2;
    rate(k) = log2(1+N^2*Pu*sum(C(k)*A(:,k).*Gamma(:,k))^2./IN(k));
end

tmp_PC = zeros(K,K);

tmp_IN = DS./(2.^rate-1);
IN_check = tmp_IN./IN;

end

