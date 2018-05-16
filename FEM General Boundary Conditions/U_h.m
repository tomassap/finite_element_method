function [ c ] = U_h(nodeFile,A,fullA,b)
% make c_g then do fullA * c_g then b-A*c_g then solve Ac_0=b-A*c_g
% then make c = c_0 + c_g and graph

% load in files
node = load(nodeFile);
[n,~]=size(fullA);
[N,~]=size(A);
c_g=zeros(n,1);
for i = 1 : n
    if node(i+1,4) == 1
       c_g(i,1) = g(node(i+1,2),node(i+1,3));
    else
       c_g(i,1) = 0; 
    end
end
Ac_g = fullA * c_g;
rhs = zeros(N,1);
count = 1;
for i = 1 : n
    if node(i+1,4) == 0
       rhs(count,1) = Ac_g(i,1);
       count = count + 1;
    end
end
rhs = b - rhs;

c_0 = A\rhs;

count2 = 1;

c = zeros(n,1);
for i = 1 : n
    if node(i+1,4) == 0
        c(i,1) = c_0(count2);
        count2 = count2 + 1;
    else
        c(i,1) = c_g(i);
    end
end

end

