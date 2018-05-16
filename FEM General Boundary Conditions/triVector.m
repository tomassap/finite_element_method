% this function finds the vector b associated with the
% respective node and ele files
% @author Andrew Tomassone

function [b] = triVector(nodeFile, eleFile, phiCo)

% load in files
node = load(nodeFile);
ele = load(eleFile);
k = 1;
h = 1;

count = 0;

% make newNode vector
newNode = zeros(node(1,1),1);
for z = 2 : node(1,1) + 1
   if node(z,4) == 0
      count = count + 1;
      newNode(z-1,1) = count;
   else
      newNode(z-1,1) = -1;
   end
end

% for each triangle
for i = 2 : ele(1,1) + 1
    A = [node(ele(i,2)+2,2),node(ele(i,2)+2,3)];
    B = [node(ele(i,3)+2,2),node(ele(i,3)+2,3)];
    C = [node(ele(i,4)+2,2),node(ele(i,4)+2,3)];
    
    a = (A(1) + B(1) - 2*C(1))/(-3*h);
    b = (A(1)-B(1))/(sqrt(3)*h);
    c = (A(2) + B(2) - 2*C(2))/(-3*h);
    d = (A(2)-B(2))/(sqrt(3)*h);
    e = (A(1) + B(1) + C(1))/3;
    ef = (A(2) + B(2) + C(2))/3;
    jac = abs(a*d-b*c);
    % for each node
    for j = 1 : 3
        % if the node is interior, run the gaussian quadrature
        if node(ele(i,j+1)+2,4) == 0 
            x1 = e; y1 = ef;
            f1 = f(x1,y1) * (phiCo(j, 1, i - 1) * x1 + phiCo(j, 2, i - 1) * y1 + phiCo(j, 3, i - 1)) * jac * 27/60;
            x2 = a * h + e; y2 = c * h + ef;
            f2 = f(x2,y2) * (phiCo(j, 1, i - 1) * x2 + phiCo(j, 2, i - 1) * y2 + phiCo(j, 3, i - 1)) * jac * 3/60;
            x3 = a * (-h/2) + b * ((h/2)*sqrt(3)) + e; y3 = c * (-h/2) + d * ((h/2) * sqrt(3)) + ef;
            f3 = f(x3,y3) * (phiCo(j, 1, i - 1) * x3 + phiCo(j, 2, i - 1) * y3 + phiCo(j, 3, i - 1)) * jac * 3/60;
            x4 = (a * (-h/2) + b * ((-h/2) * sqrt(3)) + e); y4 = c*(-h/2)+d*((-h/2)*sqrt(3))+ef;
            f4 = f(x4,y4) * (phiCo(j, 1, i - 1) * x4 + phiCo(j, 2, i - 1) * y4 + phiCo(j, 3, i - 1)) * jac * 3/60;
            x5 = a*(-h/2)+b*0+e; y5 = c*(-h/2)+d*0+ef;
            f5 = f(x5,y5) * (phiCo(j, 1, i - 1) * x5 + phiCo(j, 2, i - 1) * y5 + phiCo(j, 3, i - 1)) * jac * 8/60;
            x6 = a*(h/4)+b*((h/4)*sqrt(3))+e; y6 = c*(h/4)+d*((h/4)*sqrt(3))+ef;
            f6 = f(x6,y6) * (phiCo(j, 1, i - 1) * x6 + phiCo(j, 2, i - 1) * y6 + phiCo(j, 3, i - 1)) * jac * 8/60;
            x7 = a*(h/4)+b*((-h/4)*sqrt(3))+e; y7 = c*(h/4)+d*((-h/4)*sqrt(3))+ef;
            f7 = f(x7,y7) * (phiCo(j, 1, i - 1) * x7 + phiCo(j, 2, i - 1) * y7 + phiCo(j, 3, i - 1)) * jac * 8/60;
            s(k) = newNode(ele(i,j+1)+1);
            p(k) = 1;
            fake_b(k) = (f1+f2+f3+f4+f5+f6+f7)*(3/4)*(sqrt(3))*(h^2);
            k = k+1;
        end
    
    end
        
end
% construct the b vector using the sparse matrix system
try
    b = sparse(s,p,fake_b);
catch exception
    b = sparse(0);
end