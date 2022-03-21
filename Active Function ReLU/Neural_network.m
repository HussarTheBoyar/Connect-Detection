function [W1,W2,B1,B2] = Neural_network(Gamma,BETAA,Rate,A_max)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[M,K] = size(A_max);
class = [];
data = [];


Gamma_1 = Gamma./max(Gamma);
Rate_1  = Rate./max(Rate);
BETAA_1 = BETAA./max(BETAA);

idx = zeros(M,K);

for k = 1:K
    [~,idx(:,k)] = sort(Rate(:,k));
end

for k=1:K
    for m = 1:0.3*M
        tmp_data = [Rate_1(idx(m,k),k); Gamma_1(idx(m,k),k); BETAA_1(idx(m,k),k)];
        data = [data tmp_data];
        class = [class; A_max(idx(m,k),k)];
    end
end

% Multi layer Neural Network training
N_loop = 10;


for i = 1:length(N_loop)
    [W1,W2,B1,B2] = Multi_layer_NN(data,class,N_loop(i));
end


end

