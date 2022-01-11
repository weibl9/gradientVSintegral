function dx = glotka(t,x,A)
%LOTKA  Generailised Lotka-Volterra predator-prey model.

dx = diag([A(1,1) + A(1,2)*x(2), A(2,1) + A(2,2)*x(1)])*x; 
% dx = zeros(2,1);
% dx(1) = (A(1,1) + A(1,2)*x(1) + A(1,3)*x(2))*x(1);
% dx(2) = (A(2,1) + A(2,2)*x(1) + A(2,3)*x(2))*x(2);

end