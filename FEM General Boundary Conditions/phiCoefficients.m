% This function is used to obtain the coefficients of the
% phi functions for the respective node and ele txt files
% @author Andrew Tomassone

function [A] = phiCoefficients(nodeFile,eleFile)
nodes = load(nodeFile);
ele = load(eleFile);
A = zeros(3,3,ele(1,1));

first = [1;0;0];
second = [0;1;0];
third = [0;0;1];

for i = 2 : ele(1,1) + 1
    
    matrix = [nodes(ele(i,2)+2,2),nodes(ele(i,2)+2,3),1; nodes(ele(i,3)+2,2),nodes(ele(i,3)+2,3),1; nodes(ele(i,4)+2,2),nodes(ele(i,4)+2,3),1];
    abc1 = matrix \ first;
    abc2 = matrix \ second;
    abc3 = matrix \ third;
   
   A(1,:,i-1)=abc1;
   A(2,:,i-1)=abc2;
   A(3,:,i-1)=abc3;
    
end