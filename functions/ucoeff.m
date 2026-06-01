function [] = ucoeff()
% Purpose: To calculate the coefficients for the u equation.

% constants
global NPI NPJ Dt Cmu
% variables
global x x_u y y_v u p mu mueff SP Su F_u F_v d_u relax_u u_old rho Istart Iend ...
    Jstart Jend b aE aW aN aS aP yplus uplus k dudx dvdx

Istart = 3;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;        
        % Geometrical parameters: Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
        AREAe = AREAw;
        AREAs = x(I) - x(I-1);
        AREAn = AREAs;
        
        % eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = 0.5*(F_u(i,J)   + F_u(i-1,J))*AREAw;
        Fe = 0.5*(F_u(i+1,J) + F_u(i,J))*AREAe;
        Fs = 0.5*(F_v(I,j)   + F_v(I-1,j))*AREAs;
        Fn = 0.5*(F_v(I,j+1) + F_v(I-1,j+1))*AREAn;
        
        % eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = (mueff(I-1,J)/(x_u(i) - x_u(i-1)))*AREAw;
        De = (mueff(I,J)/(x_u(i+1) - x_u(i)))*AREAe;
        Ds = 0.25*(mueff(I-1,J) + mueff(I,J) + mueff(I-1,J-1) + mueff(I,J-1))/(y(J) - y(J-1))*AREAs;
        Dn = 0.25*(mueff(I-1,J+1) + mueff(I,J+1) + mueff(I-1,J) + mueff(I,J))/(y(J+1) - y(J))*AREAn;
        
        % The source terms
        mus = 0.25*(mueff(I-1,J) + mueff(I,J) + mueff(I-1,J-1) + mueff(I,J-1));
        mun = 0.25*(mueff(I-1,J+1) + mueff(I,J+1) + mueff(I-1,J) + mueff(I,J));
        
        if J==2 || J==NPJ+1
            if yplus(I,J) < 11.63
                SP(i,J) = -mu(I,J)*AREAs/(0.5*AREAw);
            else
                SP(i,J) = -rho(I,J) * Cmu^0.25 * k(I,J)^0.5 / uplus(I,J) *AREAs;
            end
        else
            SP(i,J) = 0.;
        end
        
        Su(i,J) = (mueff(I,J)*dudx(I,J) - mueff(I-1,J)*dudx(I-1,J)) / (x(I) - x(I-1)) + ...
            (mun*dvdx(i,j+1) - mus*dvdx(i,j)) / (y_v(j+1) - y_v(j)) - ...
            2./3. * (rho(I,J)*k(I,J) - rho(I-1,J)*k(I-1,J))/(x(I) - x(I-1));
        Su(i,J) =  Su(i,J)*AREAw*AREAs;
        
        % The coefficients (hybrid differencing scheme)
        aW(i,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(i,J) = max([-Fe, De - Fe/2, 0.]);
        if J==2
            aS(i,J) = 0.;
        else
            aS(i,J) = max([ Fs, Ds + Fs/2, 0.]);
        end
        
        if J==NPJ+1
            aN(i,J) = 0.;
        else
            aN(i,J) = max([-Fn, Dn - Fn/2, 0.]);
        end
        aPold   = 0.5*(rho(I-1,J) + rho(I,J))*AREAe*AREAn/Dt;
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(i,J) = aW(i,J) + aE(i,J) + aS(i,J) + aN(i,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % Calculation of d(i)(J) = d_u(i)(J) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32). See subroutine pccoeff.
        d_u(i,J) = AREAw*relax_u/aP(i,J);
        
        % Putting the integrated pressure gradient into the source term b(i)(J)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.       
        b(i,J) = (p(I-1,J) - p(I,J))*AREAw + Su(I,J) + aPold*u_old(i,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i)(J)       
        aP(i,J) = aP(i,J)/relax_u;
        b(i,J)  = b(i,J) + (1.0 - relax_u)*aP(i,J)*u(i,J);
        
        % now we have implemented eq. 6.36 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program.       
    end
end
end
