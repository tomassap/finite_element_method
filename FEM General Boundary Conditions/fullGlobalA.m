function [ A ] = fullGlobalA( eleFile, local )

% load in files
ele = load(eleFile);

addCount = 1;

% for each triangle
for z = 2 : ele(1,1) + 1

        i(1,addCount) = ele(z,2)+1;
        j(1,addCount) = ele(z,2)+1;
        s(1,addCount) = local(1,1,z-1);
        addCount = addCount + 1;

        i(1,addCount) = ele(z,3)+1;
        j(1,addCount) = ele(z,3)+1;
        s(1,addCount) = local(2,2,z-1); 
        addCount = addCount + 1;

        i(1,addCount) = ele(z,4)+1;
        j(1,addCount) = ele(z,4)+1;
        s(1,addCount) = local(3,3,z-1);
        addCount = addCount + 1;

        i(1,addCount) = ele(z,2)+1;
        j(1,addCount) = ele(z,3)+1;
        s(1,addCount) = local(1,2,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = ele(z,3)+1;
        j(1,addCount) = ele(z,2)+1;
        s(1,addCount) = local(2,1,z-1);
        addCount = addCount + 1;

        i(1,addCount) = ele(z,2)+1;
        j(1,addCount) = ele(z,4)+1;
        s(1,addCount) = local(1,3,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = ele(z,4)+1;
        j(1,addCount) = ele(z,2)+1;
        s(1,addCount) = local(3,1,z-1);
        addCount = addCount + 1;

        i(1,addCount) = ele(z,4)+1;
        j(1,addCount) = ele(z,3)+1;
        s(1,addCount) = local(3,2,z-1);
        addCount = addCount + 1;
        
        i(1,addCount) = ele(z,3)+1;
        j(1,addCount) = ele(z,4)+1;
        s(1,addCount) = local(2,3,z-1);
        addCount = addCount + 1;

end
% construct global A matrix with sparse matrix system
A = sparse(i,j,s);



end

% make c_g then do fullGlobalA * c_g then b-A*c_g then solve Ac_0=b-A*c_g
% then make c = c_0 + c_g and graph
