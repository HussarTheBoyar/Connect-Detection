function [rate] = SVM_rate(Mdl,Rate,Gamma,BETAA,Phi,Pu,N,C)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(Gamma);
data = [];
A_SVM = zeros(M,K);

Gamma = Gamma./max(Gamma);
Rate  = Rate./max(Rate);
BETAA = BETAA./max(BETAA);

for k=1:K
    if C(k) == 1
        for m = 1:M
            tmp_data = [Rate(m,k) Gamma(m,k) BETAA(m,k)];
            data = [data; tmp_data];
        end
    end
end

[sz1,sz2] = size(data);

if sz1*sz2 == 0
    rate = zeros(1,K);
else
    A = predict(Mdl,data);
    check_nega = find(A == -1);
    A(check_nega) = 0;

    tmp_k = 1;

    for k = 1:K
        if C(k) == 1
            A_SVM(:,k) = A( (tmp_k-1)*M+1:tmp_k*M );
            tmp_k = tmp_k+1;
        end
    end

    rate = Check_rate(Gamma,BETAA,Phi,Pu,N,A_SVM,C);
end


end

