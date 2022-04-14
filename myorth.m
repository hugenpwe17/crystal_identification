function [a1,b1,c1] = myorth(a,b,c)
a1 = a;
b1 = b - (dot(a,b)/norm(a)/norm(a))*a;
c1 = c - (dot(a1,c)/norm(a1)/norm(a1))*a1 - (dot(b1,c)/norm(b1)/norm(b1))*b1;
end