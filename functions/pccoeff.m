function [] = pccoeff()
% Purpose: To calculate the coefficients for the pc equation. 

% constants
global NPI NPJ 
% variables
global x_u y_v pc rho SP Su F_u F_v d_u d_v SMAX SAVG Istart Iend Jstart Jend ...
    b aE aW aN aS aP

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

SMAX = 0.;
SSUM = 0.;
SAVG = 0.;
	
convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces        
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.2 or fig. 6.5
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % The constant b' in eq. 6.32   
        b(I,J) = F_u(i,J)*AREAw - F_u(i+1,J)*AREAe + F_v(I,j)*AREAs - F_v(I,j+1)*AREAn;
        
        SP(I,J) = 0.;
        Su(I,J) = 0.;
        
        b(I,J) = b(I,J) + Su(I,J);
        
        SMAX = max([SMAX,abs(b(I,J))]);
        SSUM = SSUM + abs(b(I,J));
        
        % The coefficients
        aE(I,J) = (rho(I,J) + rho(I+1,J))*d_u(i+1,J)*AREAe/2;
        aW(I,J) = (rho(I-1,J) + rho(I,J))*d_u(i,J)*AREAw/2;
        aN(I,J) = (rho(I,J) + rho(I,J+1))*d_v(I,j+1)*AREAn/2;
        aS(I,J) = (rho(I,J-1) + rho(I,J))*d_v(I,j)*AREAs/2;
        
        aP(I,J) = aE(I,J) + aW(I,J) + aN(I,J) + aS(I,J) - SP(I,J);
        
        pc(I,J) = 0.;
        
        % note: At the points nearest the boundaries, some coefficients are
        % necessarily zero. For instance at I = 1 and J = 1, the coefficients
        % aS and aW are zero since they are on the outside of the calculation
        % domain. This is automatically satisfied through the initialisation
        % where d_u(i,J) and d_v(I,j) are set to zero at these points.       
    end
end
% Average error in mass balance is summed error divided by number of internal grid points
SAVG = SSUM/((Iend - Istart)*(Jend - Jstart));
end

