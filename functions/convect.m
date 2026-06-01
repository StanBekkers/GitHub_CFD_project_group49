function [] = convect()
% Purpose: To calculate the convective mass flux component pr. unit area defined in eq. 5.7 

% constants
global NPI NPJ
% variables
global x x_u y y_v u v rho F_u F_v 

for I = 2:NPI+2
    i = I;
    for J = 2:NPJ+2
        j = J;
        F_u(i,J) = (rho(I-1,J)*(x(I)-x_u(i)) + rho(I,J)*(x_u(i)-x(I-1)))*u(i,J)/(x(I)-x(I-1)); 
        F_v(I,j) = (rho(I,J-1)*(y(J)-y_v(j)) + rho(I,J)*(y_v(j)-y(J-1)))*v(I,j)/(y(J)-y(J-1)); 
    end
end
end

