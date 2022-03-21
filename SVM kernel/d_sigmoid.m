function [outputArg1] = d_sigmoid(Z)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = exp(Z)./(exp(Z)+1).^2;

end

