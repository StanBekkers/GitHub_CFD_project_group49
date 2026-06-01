function [] = calculateuplus()
% Purpose: To calculate uplus vplus and tw.

% constants
global NPI NPJ Cmu kappa ERough
% variables
global y rho k tw yplus yplus1 yplus2 uplus mu u 

viscosity();

for I = 1:NPI+1
    i = I;
    if yplus1(I,2) < 11.63
        tw(I,2)      = mu(I,2)*0.5*(u(i,2)+u(i+1,2))/(y(2)-y(1));
        yplus1(I,2)  = sqrt(rho(I,2)*abs(tw(I,2)))*(y(2)-y(1))/mu(I,2);
        yplus(I,2)   = yplus1(I,2);
        uplus(I,2)   = yplus(I,2);
    else
        tw(I,2)      = rho(I,2)*Cmu^0.25*sqrt(k(I,2))*0.5*(u(i,2)+u(i+1,2))/uplus(I,2);
        yplus1(I,2)  = sqrt(rho(I,2)*abs(tw(I,2)))*(y(2)-y(1))/mu(I,2);
        yplus(I,2)   = yplus1(I,2);
        uplus(I,2)   = log(ERough*yplus(I,2))/kappa;
    end
    
    if yplus2(I,NPJ+1) < 11.63
        tw(I,NPJ+1)      = mu(I,NPJ+1)*0.5*(u(i,NPJ+1)+u(i+1,NPJ+1))/(y(NPJ+2)-y(NPJ+1));
        yplus2(I,NPJ+1)  = sqrt(rho(I,NPJ+1)*abs(tw(I,NPJ+1)))*(y(NPJ+2)-y(NPJ+1))/mu(I,NPJ+1);
        yplus(I,NPJ+1)   = yplus2(I,NPJ+1);
        uplus(I,NPJ+1)   = yplus(I,NPJ+1);
    else
        tw(I,NPJ+1)      = rho(I,NPJ+1)*Cmu^0.25*sqrt(k(I,NPJ+1))*0.5*(u(i,NPJ+1)+u(i+1,NPJ+1))/uplus(I,NPJ+1);
        yplus2(I,NPJ+1)  = sqrt(rho(I,NPJ+1)*abs(tw(I,NPJ+1)))*(y(NPJ+2)-y(NPJ+1))/mu(I,NPJ+1);
        yplus(I,NPJ+1)   = yplus2(I,NPJ+1);
        uplus(I,NPJ+1)   = log(ERough*yplus(I,NPJ+1))/kappa;
    end          
end
