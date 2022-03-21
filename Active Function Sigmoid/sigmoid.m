function [y] = sigmoid(W,x,b)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y = 1./(1+exp(-1*(W'*x+b)));
end

