function [riemsum] = integrals(X,F)
%A function that approximates the integral of a function using
%trapezoids.
A = zeros(1,length(X)-1);
for n = 1:length(X)-1;
    A(n) = ((F(n)+F(n+1))/2)*(X(n+1)-X(n));
end
riemsum = sum(A);
end

