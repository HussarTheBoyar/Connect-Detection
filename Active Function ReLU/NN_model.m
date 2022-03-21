function [A_NN] = NN_model(W1,W2,B1,B2,data,M,K)
%UNTITLED6 Summary of this function goes here

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

