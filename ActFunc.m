% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 1 ******************************************
% This function evaluates the output of spline nonlinearity
% s is the nonlinearity input s[n]
% af is the nonlinearity structure
% x is the nonlinearity output x[n]

function [ x, af] = ActFunc( s, af)
af.s = s; % Input of the n o n l i n e a r i t y
switch ( af. aftype)

case -1 % Signed Sigmoidal
x = (2*af.Gain/(1+exp(-s*af.Slope))-af.Gain) ;
case 0 % Linear
x = s*af.Slope; % Defaut value
case 1 % Unsigned Sigmoidal
x = af.Gain/(1+exp(-s*af.Slope)) ;
case 2 % Gaussian
x = af.Gain*exp(-s^2/af.Slope ) ;
case 3 % Polynomial
x = 0 ;
for j=1:af.Pord
x = x + af.Q(j)*s.^j; % Sum of monomials
af.g(j) = s.^j;
end

case 20 % Quadratic s p l i n e
np = af.lut_len; % Number of c o n t r o l points
Su = s/af.DeltaX + ( np-1)/2; % F i r s t part of Eq . ( 7 . b )
uIndex = floor(Su) ; % Span index i in Eq . ( 7 . b )
u = Su - uIndex; % Local a b s c i s s a u in Eq . ( 7 . a )
if uIndex<1 % The index must s t a r t from 1
uIndex = 1 ;
end
if uIndex>(np-2) 
uIndex = np - 2 ; % The index cannot exceed np − 2
end
af.g = [1 u u^2]*af.C; % F i r s t part of Eq . ( 5 ) : u^T C
x = af.g*af.Q(uIndex : uIndex+2) ; % Eq . ( 5 ) : u^T C q_i
af.uIndex = uIndex; % For d e r i v a t i v e computation
af.uSpline = u; % For d e r i v a t i v e computation

otherwise % Cubic s p l i n e
np = af.lut_len; % Number of c o n t r o l points
Su = s/af.DeltaX + ( np-1)/2; % F i r s t part of Eq . ( 7 . b )
uIndex = floor(Su); % Span index i in Eq . ( 7 . b )
u = Su - uIndex; % Local a b s c i s s a u in Eq . ( 7 . a )
if uIndex<1 % The index must s t a r t from 1
uIndex = 1 ;
end
if uIndex>(np-3)  % The index cannot exceed np − 3
uIndex = np - 3;
end
af.g = [u^3 u^2 u 1]*af.C; % F i r s t part of Eq . ( 5 ) : u^T C
x = af.g*af.Q(uIndex : uIndex+3); % Eq . ( 5 ) : u^T C q_i
af.uIndex = uIndex; % For d e r i v a t i v e computation
af.uSpline = u; % For d e r i v a t i v e computation

end
end
