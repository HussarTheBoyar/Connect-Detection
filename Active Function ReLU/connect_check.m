function [get_ratio] = connect_check(Gamma,BETAA,Phi,Pu,N,C)
%M: number of AP
%K: Number of user
[M,K] = size(Gamma);

A = zeros(M,K);
A(1:1,:) = ones(1,K);
[PC2_1,IN_1,~] = PC_comp(Gamma,BETAA,Pu,Phi,N,A,C);

A = zeros(M,K);
A(2,:) = ones(1,K);
[PC2_2,IN_2,~] = PC_comp(Gamma,BETAA,Pu,Phi,N,A,C);

A = zeros(M,K);
A(1:2,:) = ones(2,K);
[PC2,IN,~] = PC_comp(Gamma,BETAA,Pu,Phi,N,A,C);

PC_check = PC2-PC2_1-PC2_2;
tmp_PC = PC2_1.*PC2_2;

for k=1:K
    tmp_PC(k,k) = 0;
end

IN_check = (IN_1+IN_2+2*N^2*Pu*sum(tmp_PC))./IN;

PC_comp(Gamma,BETAA,Pu,Phi,N,A,C);


end

