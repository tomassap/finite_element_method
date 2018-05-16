function [ s ] = Heron( A,B,C )
a = sqrt((A(2) - B(2))^2 + (A(1) - B(1))^2);
b = sqrt((B(2) - C(2))^2 + (B(1) - C(1))^2);
c = sqrt((A(2) - C(2))^2 + (A(1) - C(1))^2);
s = (a+b+c)/2;
s = sqrt(s*(s-a)*(s-b)*(s-c));

end

