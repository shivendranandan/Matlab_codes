% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 6 ******************************************

% Nonlinear Function Implementation


function value = G1_FUNC( X, G, S, ty)

switch ty
    case -1
        value = 2*G/(1+exp(-X*S)) - G; % Signed Sigmoid
    case 0
        value = X*S; % Linear
    case 1
        value = G/(1+exp(-X*S)); % Unsigned Sigmoidal
    case 2
        value = G* exp(-(X*X)/5.0) ; % Gaussian
    case 3
        value = 2*G*(rand - 0.5); % Random
    otherwise
        value = 0 ;
        
end
end % End function FUNC( ) ===================================================