function [] = ucoeff()
% Purpose: To calculate the coefficients for the u equation.

% Constants
global NPI NPJ LARGE XMAX YMAX
% variables
global x x_u y y_v u p mueff SP Su F_u F_v d_u relax_u Istart Iend Jstart Jend ...
      b aE aW aN aS aP
global h_base_frac l_base_frac
 
h_base_frac = 2/10;
l_base_frac = 3/10;

Istart = 3;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

convect();
layout_wall = Walls(Istart, Iend, Jstart, Jend, NPI, NPJ, h_base_frac);
layout_fins = TriangleFin(Istart, Iend, Jstart, Jend, NPI, NPJ, l_base_frac, h_base_frac);
cooler_layout = layout_wall | layout_fins;   % or use + and threshold if needed

% cell_area = (XMAX*YMAX)/(NPI*NPJ);
% num_ones = sum(cooler_layout(:));
% fin_area = num_ones* cell_area
% 
% figure
% subplot(1,2,1)
% imagesc(layout_wall', [0 1])
% colormap(gray)
% axis equal tight
% colorbar
% title('layout\_wall (j,i)')
% xlabel('j'); ylabel('i')
% 
% subplot(1,2,2)
% imagesc(layout_fins', [0 1])
% colormap(gray)
% axis equal tight
% colorbar
% title('layout\_fins (j,i)')
% xlabel('j'); ylabel('i')
% 
% figure
% imagesc(cooler_layout', [0 1])
% colormap(gray)
% axis equal tight
% colorbar
% title('combined layout (j,i)')
% xlabel('j'); ylabel('i')

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
        Fw = ((F_u(i,J)   + F_u(i-1,J))/2)*AREAw;
        Fe = ((F_u(i+1,J) + F_u(i,J))/2)*AREAe;
        Fs = ((F_v(I,j)   + F_v(I-1,j))/2)*AREAs;
        Fn = ((F_v(I,j+1) + F_v(I-1,j+1))/2)*AREAn;
        
        % eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = (mueff(I-1,J)/(x_u(i) - x_u(i-1)))*AREAw;
        De = (mueff(I,J)/(x_u(i+1) - x_u(i)))*AREAe;
        Ds = ((mueff(I-1,J) + mueff(I,J) + mueff(I-1,J-1) + mueff(I,J-1))/(4*(y(J) - y(J-1))))*AREAs;
        Dn = ((mueff(I-1,J+1) + mueff(I,J+1) + mueff(I-1,J) + mueff(I,J))/(4*(y(J+1) - y(J))))*AREAn;
        
        % The source terms
        SP(i,J) = 0.;
        Su(i,J) = 0.;
        
         % u can be fixed to zero by setting SP to a very large value    
        if (cooler_layout(i,j) == 1)
            SP(i,J) = -LARGE;
        end

        % The coefficients (hybrid differencing scheme)
        aW(i,J) = max([ Fw, Dw + Fw/2, 0.]);
        aE(i,J) = max([-Fe, De - Fe/2, 0.]);
        aS(i,J) = max([ Fs, Ds + Fs/2, 0.]);
        aN(i,J) = max([-Fn, Dn - Fn/2, 0.]);
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(i,J) = aW(i,J) + aE(i,J) + aS(i,J) + aN(i,J) + Fe - Fw + Fn - Fs - SP(I,J);
        
        % Calculation of d(i)(J) = d_u(i)(J) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32). See subroutine pccoeff.
        d_u(i,J) = AREAw*relax_u/aP(i,J);
        
        % Putting the integrated pressure gradient into the source term b(i)(J)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.     
        b(i,J) = (p(I-1,J) - p(I,J))*AREAw + Su(I,J);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i)(J)
        aP(i,J) = aP(i,J) / relax_u;
        b (i,J) = b(i,J) + (1 - relax_u)*aP(i,J)*u(i,J);
        
        % now we have implemented eq. 6.36 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program. 
    end
end



end

