function fi = solve(fi,b,aE,aW,aN,aS,aP)
% Purpose: To solve the algebraic equation 7.7.

% variables
global Istart Iend Jstart Jend 

% TDMA along a horizontal row from west to east. equation to solve:
% - aW*fiW + aP*fiP - aE*fiE = aS*fiS + aN*fiN + b

% equivalences with variables in eq. 7.1-7.6:
% BETA = aW(I,J) Def. in eq. 7.2
% D    = aP(I,J) Def. in eq. 7.2
% ALFA = aE(I,J) Def. in eq. 7.2
% A    = Ari(I)	 Def. in eq. 7.6b
% C    = Cri	 The right side aSAVGed temporarily known (see eq. 7.8)
% C'   = Cmri(I)  Def. in eq. 7.6c
% b    = b(I,J)	 Def. in eq. 7.7

space = max([(Iend - Istart + 3),(Jend - Jstart + 3)]);
Ari   = zeros(1,space);
Cmri  = zeros(1,space);

% Solving the (e-w) lines from the south
for J = Jstart:Jend
    % At the inlet boundary:
    Ari (Istart-1) = 0;
    Cmri(Istart-1) = fi(Istart-1,J);
    
    for I = Istart:Iend % Forward substitution
        Ari(I)  = aE(I,J)/(aP(I,J) - aW(I,J)*Ari(I-1)); % eq. 7.6b
        Cri     = aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J);
        Cmri(I) = (aW(I,J)*Cmri(I-1) + Cri)/(aP(I,J) - aW(I,J)*Ari(I-1)); % eq. 7.6c
    end
    
    for I = Iend:-1:Istart  % Back substitution
        fi(I,J) = Ari(I)*fi(I+1,J) + Cmri(I); % eq. 7.6a
    end
end

% Solving the (e-w) lines from the north
for J = Jend-1:-1:Jstart
    % At the inlet boundary:
    Ari (Istart-1) = 0;
    Cmri(Istart-1) = fi(Istart-1,J);
    
    for I = Istart:Iend % Forward substitution
        Ari(I)  = aE(I,J)/(aP(I,J) - aW(I,J)*Ari(I-1)); % eq. 7.6b
        Cri     = aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J);
        Cmri(I) = (aW(I,J)*Cmri(I-1) + Cri)/(aP(I,J) - aW(I,J)*Ari(I-1));  % eq. 7.6c
    end
    
    for I = Iend:-1:Istart  % Back substitution
        fi(I,J) = Ari(I)*fi(I+1,J) + Cmri(I); % eq. 7.6a
    end
end

% TDMA along a vertical column from south to north. equation to solve:
% - aS*fiW + aP*fiP - aN*fiE = aW*fiS + aE*fiN + b (eq. 7.8)

% equivalences with variables in eq. 7.1-7.6:
% BETA = aS(I,J) Def. in eq. 7.2
% D    = aP(I,J) Def. in eq. 7.2
% ALFA = aN(I,J) Def. in eq. 7.2
% A    = Ari(I)	 Def. in eq. 7.6b
% C    = Cri      The right side aSAVGed temporarily known (see eq. 7.8)
% C'  = Cmri(I)  Def. in eq. 7.6c
% b    = b(I,J)	 Def. in eq. 7.7

% Solving (n-s) lines from the west
for I = Istart:Iend
    % At the bottom boundary:
    Ari(Jstart-1) = 0;
    Cmri(Jstart-1) = fi(I,Jstart-1);
    
    for J = Jstart:Jend % Forward substitution
        Ari(J)  = aN(I,J)/(aP(I,J) - aS(I,J)*Ari(J-1)); % eq. 7.6b
        Cri     = aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) + b(I,J);
        Cmri(J) = (aS(I,J)*Cmri(J-1) + Cri)/(aP(I,J) - aS(I,J)*Ari(J-1)); % eq. 7.6c
    end
    
    for J = Jend:-1:Jstart % Back substitution
        fi(I,J) = Ari(J)*fi(I,J+1) + Cmri(J); % eq. 7.6a
    end
end

% Solving (n-s) lines from the east
for I = Iend-1:-1:Istart
    % At the bottom boundary:
    Ari(Jstart-1) = 0;
    Cmri(Jstart-1) = fi(I,Jstart-1);
    
    for J = Jstart:Jend % Forward substitution
        Ari(J)  = aN(I,J)/(aP(I,J) - aS(I,J)*Ari(J-1)); % eq. 7.6b
        Cri     = aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) + b(I,J);
        Cmri(J) = (aS(I,J)*Cmri(J-1) + Cri)/(aP(I,J) - aS(I,J)*Ari(J-1)); % eq. 7.6c
    end
    
    for J = Jend:-1:Jstart % Back substitution
        fi(I,J) = Ari(J)*fi(I,J+1) + Cmri(J); % eq. 7.6a
    end
end
end
