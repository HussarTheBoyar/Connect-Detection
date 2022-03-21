function [A] = ReLU(W,x,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
A = max(0,W'*x+b);
end

