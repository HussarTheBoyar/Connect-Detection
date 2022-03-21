function [A_NN] = NN_model_new(W1,W2,B1,B2,M,K,Gamma,BETAA,Rate,C)
%UNTITLED6 Summary of this function goes here

Gamma = Gamma./max(Gamma);
Rate  = Rate./max(Rate);
BETAA = BETAA./max(BETAA);

data = [];

for k=1:K
    if C(k) == 1
            for m = 1:M
                tmp_data = [Rate(m,k); Gamma(m,k); BETAA(m,k)];
                data = [data tmp_data];
            end
    end
end

%Hidden layer state
A1 = ReLU(W1,data,B1);
%Output Layer
Y_output = sigmoid(W2,A1,B2);

threshold = mean(Y_output);
Y_hat = Y_output>=threshold;

A_NN = zeros(M,K);

for k=1:K
    A_NN(:,k) = Y_hat( (k-1)*M+1:k*M );
end

end

