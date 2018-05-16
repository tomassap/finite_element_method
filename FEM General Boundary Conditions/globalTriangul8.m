% this function uses the local A matrix in the triangul8 function to
% find the global A matrix
% @author Andrew Tomassone

function [global_matrix] = globalTriangul8(nodeFile, eleFile, local)

% load in files
node = load(nodeFile);
ele = load(eleFile);
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
addCount = 1;

% for each triangle
for z = 2 : ele(1,1) + 1
    if newNode(ele(z,2)+1) > 0
        i(1,addCount) = newNode(ele(z,2)+1);
        j(1,addCount) = newNode(ele(z,2)+1);
        s(1,addCount) = local(1,1,z-1);
        addCount = addCount + 1;
    end
    if newNode(ele(z,3)+1) > 0
        i(1,addCount) = newNode(ele(z,3)+1);
        j(1,addCount) = newNode(ele(z,3)+1);
        s(1,addCount) = local(2,2,z-1); 
        addCount = addCount + 1;
    end
    if newNode(ele(z,4)+1) > 0 
        i(1,addCount) = newNode(ele(z,4)+1);
        j(1,addCount) = newNode(ele(z,4)+1);
        s(1,addCount) = local(3,3,z-1);
        addCount = addCount + 1;
    end
    if newNode(ele(z,2)+1) > 0 && newNode(ele(z,3)+1) > 0
        i(1,addCount) = newNode(ele(z,2)+1);
        j(1,addCount) = newNode(ele(z,3)+1);
        s(1,addCount) = local(1,2,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = newNode(ele(z,3)+1);
        j(1,addCount) = newNode(ele(z,2)+1);
        s(1,addCount) = local(2,1,z-1);
        addCount = addCount + 1;
    end
    if newNode(ele(z,2)+1) > 0 && newNode(ele(z,4)+1) > 0
        i(1,addCount) = newNode(ele(z,2)+1);
        j(1,addCount) = newNode(ele(z,4)+1);
        s(1,addCount) = local(1,3,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = newNode(ele(z,4)+1);
        j(1,addCount) = newNode(ele(z,2)+1);
        s(1,addCount) = local(3,1,z-1);
        addCount = addCount + 1;
    end
    if newNode(ele(z,3)+1) > 0 && newNode(ele(z,4)+1) > 0
        i(1,addCount) = newNode(ele(z,4)+1);
        j(1,addCount) = newNode(ele(z,3)+1);
        s(1,addCount) = local(3,2,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = newNode(ele(z,3)+1);
        j(1,addCount) = newNode(ele(z,4)+1);
        s(1,addCount) = local(2,3,z-1);
        addCount = addCount + 1;
    end
end
% construct global A matrix with sparse matrix system
try
    global_matrix = sparse(i,j,s);
catch exception
    global_matrix = sparse(0);
end
    