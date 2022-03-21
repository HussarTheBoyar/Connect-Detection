function [W1,W2,B1,B2] = Multi_layer_NN(data,class,N_loop)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Innitial
X = data;
[row,~] = size(X);
N = 5000;
W1 = rand(row,N);
W2 = rand(N,1);
B1 = rand(N,1);
B2 = rand;
Y = class;
alpha = 0.0005;

A_1 = ReLU(W1,X(:,1),B1);
Z_1 = W1'*X(:,1)+B1;
Y_hat = sigmoid(W2,A_1,B2);

%Processing

for count = 1:N_loop
    for i = 1:length(X)
%Update parameter        
        E2 = abs(Y_hat - Y(i));
        delta_W2 = A_1*E2';
        delta_B2 = sum(E2);
        
        W2 = W2-alpha*delta_W2;
        B2 = B2-alpha*delta_B2;
        
        E1 = (W2*E2).*d_ReLU(Z_1);
        delta_W1 = X(:,i)*E1';
        delta_B1 = sum(E1);
        
        W1 = W1-alpha*delta_W1;
        B1 = B1-alpha*delta_B1;
%Update state
        
        A_1 = ReLU(W1,X(:,1),B1);
        Z_1 = W1'*X(:,1)+B1;
        Y_hat = sigmoid(W2,A_1,B2);
        
    end

    swap = randperm(length(X));
    tmp_data = X(:,swap);
    tmp_Y = Y(swap);
    X = tmp_data;
    Y = tmp_Y;
    
end

end

