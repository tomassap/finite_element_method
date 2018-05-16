
% This script calculates L_2 error and more corresponding to the
% associated finite element method programs.
%
%
% @author Andrew Tomassone


for i = 1:7
    nodefile = load(sprintf('Lmeshnode%d.txt', i));
    elefile = load(sprintf('Lmeshele%d.txt', i));
    phiCo = phiCoefficients(sprintf('Lmeshnode%d.txt', i), sprintf('Lmeshele%d.txt', i));
    count = 0;
    totalError = 0;
    totalError_grad_u = 0;
    totalError_v_h = 0;
    totalError_grad_v = 0;
    c_g = zeros(nodefile(1, 1));
    % make c_g
    for j = 2: nodefile(1, 1) + 1
        if nodefile(j, 4) == 1
            c_g(j - 1) = g(nodefile(j, 2), nodefile(j, 3));
        else 
            c_g(j - 1) = 0;
        end
    end
    
    
    newA = fullGlobalA(sprintf('Lmeshele%d.txt', i), triangul8(sprintf('Lmeshnode%d.txt', i),sprintf('Lmeshele%d.txt', i))); 
    A = globalTriangul8(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i),triangul8(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i)));
    b = triVector(sprintf('Lmeshnode%d.txt',i),sprintf('Lmeshele%d.txt',i), phiCo);
     u_h = U_h(sprintf('Lmeshnode%d.txt',i),A,newA,b);
    % x = A \ b;
    product = newA * c_g; % must shorten this
    short_c_g = zeros(size(b));
    count = 1;
    for j = 1: nodefile(1, 1)
        if nodefile(j + 1, 4) == 0
            short_c_g(count) = product(j);
            count = count + 1;
        end
    end
    
    rhs = b - short_c_g;
    c_0 = A \ rhs;
    
    % add c_0 and c_g
    count = 1;
    for h = 1: nodefile(1, 1)
        if nodefile(h + 1, 4) == 0
            c_g(h) = c_0(count);
            count = count + 1;
        end
    end
    x = c_g;
    
    for z = 2: nodefile(1,1) + 1
        v(z - 1) = u(nodefile(z, 2), nodefile(z, 3));
    end 

    % loop through the triangle
    for j = 2: elefile(1, 1) + 1
        % do quadrature
        A(1) = nodefile(elefile(j, 2) + 2, 2);
        A(2) = nodefile(elefile(j, 2) + 2, 3);
 
        B(1) = nodefile(elefile(j, 3) + 2, 2);
        B(2) = nodefile(elefile(j, 3) + 2, 3);
 
        C(1) = nodefile(elefile(j, 4) + 2, 2);
        C(2) = nodefile(elefile(j, 4) + 2, 3);
 
        h = 10^14;
 
        a = (A(1) + B(1) - 2*C(1))/(-3*h);
        b = (A(1) - B(1))/(sqrt(3) * h);
        c = (A(2) + B(2) - 2*C(2))/(-3*h);
        d = (A(2) - B(2))/(sqrt(3) * h);
        e = (A(1) + B(1) + C(1))/3;
        eff = (A(2) + B(2) + C(2))/3; % avoid name clash with f.m
 
        jacob = abs(a * d - b * c);
 
        % initalize f and x's and y's
        f1 = 0;
        f2 = 0;
        f3 = 0;
        f4 = 0;
        f5 = 0;
        f6 = 0;
        f7 = 0;
        f1_2 = 0;
        f2_2 = 0;
        f3_2 = 0;
        f4_2 = 0;
        f5_2 = 0;
        f6_2 = 0;
        f7_2 = 0;
        u_h_x = 0;
        u_h_y = 0;
        v_h_x = 0;
        v_h_y = 0;
        x1 = e;
        y1 = eff;
        x2 = a*h+e; 
        y2 = c*h+eff;
        x3 = a * (-h/2) + b * (h * sqrt(3) / 2) + e; 
        y3 = c * (-h/2) + d * (h * sqrt(3) / 2) + eff;
        x4 = a * (-h/2) + b * (-h * sqrt(3) / 2) + e;
        y4 = c * (-h/2) + d * (-h * sqrt(3) / 2) + eff;
        x5 = a * (-h/2) + b * 0 + e;
        y5 = c * (-h/2) + d * 0 + eff;
        x6 = a * (h/4) + b * ((h / 4) * sqrt(3)) + e;
        y6 = c * (h/4) + d * ((h / 4) * sqrt(3)) + eff;
        x7 = a * (h/4) + b * (-h / 4 * sqrt(3)) + e;
        y7 = c * (h/4) + d * (-h / 4 * sqrt(3)) + eff;
        % calculate uh
        for k = 1: 3
                % make derivitives for each gradient
                u_h_x = u_h_x + phiCo(k, 1, j - 1) * u_h(elefile(j, k + 1) + 1);
                u_h_y = u_h_y + phiCo(k, 2, j - 1) * u_h(elefile(j, k + 1) + 1);
                v_h_x = v_h_x + phiCo(k, 1, j - 1) * v(elefile(j, k + 1) + 1);
                v_h_y = v_h_y + phiCo(k, 2, j - 1) * v(elefile(j, k + 1) + 1);
                % make u_h's
                f1 = f1 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x1 + phiCo(k, 2, j - 1) * y1 + phiCo(k, 3, j - 1));
                f2 = f2 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x2 + phiCo(k, 2, j - 1) * y2 + phiCo(k, 3, j - 1));
                f3 = f3 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x3 + phiCo(k, 2, j - 1) * y3 + phiCo(k, 3, j - 1));
                f4 = f4 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x4 + phiCo(k, 2, j - 1) * y4 + phiCo(k, 3, j - 1));
                f5 = f5 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x5 + phiCo(k, 2, j - 1) * y5 + phiCo(k, 3, j - 1));
                f6 = f6 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x6 + phiCo(k, 2, j - 1) * y6 + phiCo(k, 3, j - 1));
                f7 = f7 + x(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x7 + phiCo(k, 2, j - 1) * y7 + phiCo(k, 3, j - 1));
                % make v_h's
                f1_2 = f1_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x1 + phiCo(k, 2, j - 1) * y1 + phiCo(k, 3, j - 1));
                f2_2 = f2_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x2 + phiCo(k, 2, j - 1) * y2 + phiCo(k, 3, j - 1));
                f3_2 = f3_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x3 + phiCo(k, 2, j - 1) * y3 + phiCo(k, 3, j - 1));
                f4_2 = f4_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x4 + phiCo(k, 2, j - 1) * y4 + phiCo(k, 3, j - 1));
                f5_2 = f5_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x5 + phiCo(k, 2, j - 1) * y5 + phiCo(k, 3, j - 1));
                f6_2 = f6_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x6 + phiCo(k, 2, j - 1) * y6 + phiCo(k, 3, j - 1));
                f7_2 = f7_2 + v(elefile(j, k + 1) + 1) * (phiCo(k, 1, j - 1) * x7 + phiCo(k, 2, j - 1) * y7 + phiCo(k, 3, j - 1));
        end
        % do (u - uh)^2, multiply jacobian and weight
         f1 = (u(x1, y1) - f1)^2 * jacob * 27/60;
         f2 = (u(x2, y2) - f2)^2 * jacob * 3/60;
         f3 = (u(x3, y3) - f3)^2 * jacob * 3/60;
         f4 = (u(x4, y4) - f4)^2 * jacob * 3/60;
         f5 = (u(x5, y5) - f5)^2 * jacob * 8/60;
         f6 = (u(x6, y6) - f6)^2 * jacob * 8/60;
         f7 = (u(x7, y7) - f7)^2 * jacob * 8/60;
        
        % do (grad(u) - grad(u_h))^2, multiply jacobian and weight
        % find u_x and u_y by taking those derivitives of u and using them
        % as functions
         f1_grad_u = (grad_u(x1, y1) - [u_h_x, u_h_y]).^2 * jacob * 27/60;
         f2_grad_u = (grad_u(x2, y2) - [u_h_x, u_h_y]).^2 * jacob * 3/60;
         f3_grad_u = (grad_u(x3, y3) - [u_h_x, u_h_y]).^2 * jacob * 3/60;
         f4_grad_u = (grad_u(x4, y4) - [u_h_x, u_h_y]).^2 * jacob * 3/60;
         f5_grad_u = (grad_u(x5, y5) - [u_h_x, u_h_y]).^2 * jacob * 8/60;
         f6_grad_u = (grad_u(x6, y6) - [u_h_x, u_h_y]).^2 * jacob * 8/60;
         f7_grad_u = (grad_u(x7, y7) - [u_h_x, u_h_y]).^2 * jacob * 8/60;

         % do (u - vh)^2, multiply jacobian and weight
         f1_2 = (u(x1, y1) - f1_2)^2 * jacob * 27/60;
         f2_2 = (u(x2, y2) - f2_2)^2 * jacob * 3/60;
         f3_2 = (u(x3, y3) - f3_2)^2 * jacob * 3/60;
         f4_2 = (u(x4, y4) - f4_2)^2 * jacob * 3/60;
         f5_2 = (u(x5, y5) - f5_2)^2 * jacob * 8/60;
         f6_2 = (u(x6, y6) - f6_2)^2 * jacob * 8/60;
         f7_2 = (u(x7, y7) - f7_2)^2 * jacob * 8/60;

        % do (grad(u) - grad(u_h))^2, multiply jacobian and weight
        % find u_x and u_y by taking those derivitives of u and using them
        % as functions
         f1_grad_v = (grad_u(x1, y1) - [v_h_x, v_h_y]).^2 * jacob * 27/60;
         f2_grad_v = (grad_u(x2, y2) - [v_h_x, v_h_y]).^2 * jacob * 3/60;
         f3_grad_v = (grad_u(x3, y3) - [v_h_x, v_h_y]).^2 * jacob * 3/60;
         f4_grad_v = (grad_u(x4, y4) - [v_h_x, v_h_y]).^2 * jacob * 3/60;
         f5_grad_v = (grad_u(x5, y5) - [v_h_x, v_h_y]).^2 * jacob * 8/60;
         f6_grad_v = (grad_u(x6, y6) - [v_h_x, v_h_y]).^2 * jacob * 8/60;
         f7_grad_v = (grad_u(x7, y7) - [v_h_x, v_h_y]).^2 * jacob * 8/60;
 
        % add to total error of triangles
        totalError = totalError + (3/4) * (h^2 * sqrt(3)) * (f1 + f2 + f3 + f4 + f5 + f6 + f7);
        totalError_grad_u = totalError_grad_u + (3/4) * (h^2 * sqrt(3)) * (f1_grad_u + f2_grad_u + f3_grad_u + f4_grad_u + f5_grad_u + f6_grad_u + f7_grad_u);
        totalError_v_h = totalError_v_h + (3/4) * (h^2 * sqrt(3)) * (f1_2 + f2_2 + f3_2 + f4_2 + f5_2 + f6_2 + f7_2);
        totalError_grad_v = totalError_grad_u + (3/4) * (h^2 * sqrt(3)) * (f1_grad_v + f2_grad_v + f3_grad_v + f4_grad_v + f5_grad_v + f6_grad_v + f7_grad_v);
    
    end
   fprintf('Approximation %d error ||u - u_h||: %e\n', i, sqrt(totalError));
   fprintf('Approximation %d error ||grad(u) - grad(u_h)||: %e\n', i, sqrt(totalError_grad_u(1, 1)));
   fprintf('Approximation %d error ||u - v_h||: %e\n', i, sqrt(totalError));
   fprintf('Approximation %d error ||grad(u) - grad(v_h)||: %e\n', i, sqrt(totalError_grad_v(1, 1)));
end