% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 7 ******************************************

% For Calculation of Table Length


function Table_Length = G1_TabFuncLen( DX, G, S, ty)
ii = 0 ;
X = 0 ;
if (ty~=2 && ty~=3)
    F = 0.0 ;
    crtGain = G - 0.005*G; % Max a t c func value −−−−−−−−−−−−−−−−−−−−−
    while ( F<crtGain)
        F = G1_FUNC( X, G, S, ty) ;
        X = X + DX;
        ii = ii + 1 ;
    end
elseif (ty==3)
    ii = 11 ;
else
    crtGain = 0.005*G; % Gaussian
    F = G;
    while ( F>crtGain )
        ii = ii + 1 ;
        F = FUNC( X, G, S, ty) ;
        X = X + DX;
    end
end
Table_Length = ii*2 + 1; % always odd
end

% End function TabFuncLen ( ) =============================================