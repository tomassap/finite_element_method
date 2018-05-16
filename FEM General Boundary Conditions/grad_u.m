function grad = grad_u(x,y)
grad = [y*(2*x-1)*(y-1), x*(x-1)*(2*y-1)];
end