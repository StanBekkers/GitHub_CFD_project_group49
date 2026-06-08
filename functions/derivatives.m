function [] = derivatives()
% Purpose: To calculate derivatives
% Fix: E(i,j) is on the staggered velocity grid, indexed up to NPI+2, NPJ+2.
% E2(I,J) averages 4 E values at corners — the loop must stop at NPI and NPJ
% to avoid accessing E(NPI+2, NPJ+2) which is out of bounds.

% global variables
global x x_u y y_v u v dudx dudy dvdx dvdy E E2
% global constants
global NPI NPJ

% Initialise to zero so boundary cells are clean
dudx = zeros(NPI+2, NPJ+2);
dudy = zeros(NPI+2, NPJ+2);
dvdx = zeros(NPI+2, NPJ+2);
dvdy = zeros(NPI+2, NPJ+2);
E    = zeros(NPI+2, NPJ+2);
E2   = zeros(NPI+2, NPJ+2);

% First pass: compute strain-rate components on staggered grid
for I = 2:NPI+1
    i = I;
    for J = 2:NPJ+1
        j = J;
        dudx(I,J) = (u(i+1,J) - u(i,J))   / (x_u(i+1) - x_u(i));
        dudy(i,j) = (u(i,J)   - u(i,J-1)) / (y(J)     - y(J-1));
        dvdx(i,j) = (v(I,j)   - v(I-1,j)) / (x(I)     - x(I-1));
        dvdy(I,J) = (v(I,j+1) - v(I,j))   / (y_v(j+1) - y_v(j));
        E(i,j)    = (dudy(i,j))^2 + (dvdx(i,j))^2 + 2*dudy(i,j)*dvdx(i,j);
    end
end

% Second pass: E2 needs E at 4 corners (i,j), (i+1,j), (i,j+1), (i+1,j+1)
% So i and j can only go up to NPI+1 and NPJ+1 safely.
% Loop over I=2:NPI, J=2:NPJ to keep corner indices within bounds.
for I = 2:NPI
    i = I;
    for J = 2:NPJ
        j = J;
        E2(I,J) = dudx(I,J)^2 + dvdy(I,J)^2 + ...
                  0.25*(E(i,j) + E(i+1,j) + E(i,j+1) + E(i+1,j+1));
    end
end
% Fill edge column and row by copying nearest interior value
E2(NPI+1, 2:NPJ)   = E2(NPI,   2:NPJ);
E2(2:NPI,  NPJ+1)  = E2(2:NPI,  NPJ);
E2(NPI+1,  NPJ+1)  = E2(NPI,    NPJ);

end