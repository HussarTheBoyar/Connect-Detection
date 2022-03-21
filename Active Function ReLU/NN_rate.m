function [rate,A_NN] = NN_rate(W1,W2,B1,B2,Rate,Gamma,BETAA,Phi,Pu,N,C)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(Gamma);

data = [];

Gamma_1 = Gamma./max(Gamma);
Rate_1  = Rate./max(Rate);
BETAA_1 = BETAA./max(BETAA);

for k=1:K
    for m = 1:M
        tmp_data = [Rate_1(m,k); Gamma_1(m,k); BETAA_1(m,k)];
        data = [data tmp_data];
    end
end
[sz1,sz2] = size(data);

if sz1*sz2 == 0
    rate = zeros(1,K);
else
    [A_NN] = NN_model(W1,W2,B1,B2,data,M,K);

    rate = Check_rate(Gamma,BETAA,Phi,Pu,N,A_NN,C);
end


end

