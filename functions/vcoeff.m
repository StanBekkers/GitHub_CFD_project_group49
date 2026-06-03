function [] = vcoeff()
% Purpose: To calculate the coefficients for the v equation.

% constants
global NPI NPJ 
% variables
global x x_u y y_v v p mu SP Su F_u F_v d_v relax_v Istart Iend Jstart Jend ...
    b aE aW aN aS aP

Istart = 2;
Iend = NPI+1;
Jstart = 3;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y(J) - y(J-1); % See fig. 6.4
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;
        
        % eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = ((F_u(i,J)   + F_u(i,J-1))/2)*AREAw;
        Fe = ((F_u(i+1,J) + F_u(i+1,J-1))/2)*AREAe;
        Fs = ((F_v(I,j)   + F_v(I,j-1))/2)*AREAs;
        Fn = ((F_v(I,j)   + F_v(I,j+1))/2)*AREAn;
        
        % eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = ((mu(I-1,J-1) + mu(I,J-1) + mu(I-1,J) + mu(I,J))/(4*(x(I) - x(I-1))))*AREAw;
        De = ((mu(I,J-1) + mu(I+1,J-1) + mu(I,J) + mu(I+1,J))/(4*(x(I+1) - x(I))))*AREAe;
        Ds =  (mu(I,J-1)/(y_v(j) - y_v(j-1)))*AREAs;
        Dn =  (mu(I,J)/(y_v(j+1) - y_v(j)))*AREAn;
        
        % The source terms
        SP(I,j) = 0.;
        Su(I,j) = 0.;
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);
        
        % transport of v through the baffles can be switched off by setting the coefficients to zero
        
        %lower walls: 
         if (J < ceil((NPJ+1)/6)) 
            aW(I,j) = 0;
            aE(I,j) = 0;
         end
        %upper walls:
        if (J > ceil(5*(NPJ+1)/6)) 
            aW(I,j) = 0;
            aS(I,j) = 0; 
        end

 N_t_points_u = ceil(5*(NPI+1)/10) - ceil(4*(NPI+1)/10);
     
        %first triangles: 
        if (ceil(4*(NPI+1)/10)-1<= i && i <=(ceil(5*(NPI+1)/10)) && J < ceil( (1/ N_t_points_u) * ceil(3*(NPJ+1)/10) + ceil((NPJ+1)/6)))
            aE(I,j) = 0;
        end
        if (ceil(4*(NPI+1)/10)<= i && i <=(ceil(5*(NPI+1)/10)) && J < ceil( (1/ N_t_points_u) * ceil(3*(NPJ+1)/10) + ceil((NPJ+1)/6)))
            aW(I,j) = 0;
        end
        if (ceil(4*(NPI+1)/10)-1<= i && i <=(ceil(5*(NPI+1)/10)) && J > ceil((-1/N_t_points_u)*ceil(7*(NPJ+1)/10) +ceil(5*(NPJ+1)/6)))
            aE(I,j) = 0;
        end
        if (ceil(4*(NPI+1)/10)<= i && i <=(ceil(5*(NPI+1)/10)) && J > ceil((-1/N_t_points_u)*ceil(7*(NPJ+1)/10) +ceil(5*(NPJ+1)/6)))
            aW(I,j) = 0;
        end
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,j) = aW(I,j) + aE(I,j) + aS(I,j) + aN(I,j) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Calculation of d(I,j) = d_v(I,j) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32) (see subroutine pccoeff).
        d_v(I,j) = AREAs*relax_v/aP(I,j);
        
        % Putting the integrated pressure gradient into the source term b(I,j)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.
        b(I,j) = (p(I,J-1) - p(I,J))*AREAs + Su(I,j);
        
        % Introducing relaxation by eq. 6.37 . and putting also the last
        % term on the right side into the source term b(i,J)
        aP(I,j) = aP(I,j)/relax_v;
        b(I,j)  = b(I,j) + (1 - relax_v)*aP(I,j)*v(I,j);
        
        % now we have implemented eq. 6.37 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program.        
    end
end
end

