% this script produces graphs for the triangulations for the first 6 level meshes
% run the script to output the figures, which will initialize to a top
% level view - to rotate, press the counter-clockwise blue arrow
% @author Andrew Tomassone

% update i to iterate from 1 : 8 if your cpu can handle larger meshes
n = 6;
for i = 1 : n
    % load in files, obtain A and b then solve for x
    ele = load(sprintf('Lmeshele%d.txt',i));
    node = load(sprintf('Lmeshnode%d.txt',i));
    A = globalTriangul8(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i),triangul8(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i)));
    phis = phiCoefficients(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i));
    b = triVector(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i), phis);
    fullA = fullGlobalA(sprintf('Lmeshele%d.txt',i),triangul8(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i)));
    u_h = U_h(sprintf('Lmeshnode%d.txt',i),A,fullA,b);
    
    % for each triangle
    for j = 2 : ele(1,1) + 1
        
        % for each node
        for k = 1 : 3
            xes(k) = node(ele(j,k+1)+2,2);
            yes(k) = node(ele(j,k+1)+2,3);
            zes(k) = u_h(ele(j,k+1)+1,1); 
        end
        
        % use patch to graph approximation and hold as you iterate through
        % the triangles
        patch(xes,yes,zes,[.9,1,1]);
        hold on;  
    end
    
    % if this is not the last mesh, close the figure
    if (i < n)
        figure();
    end

    
end