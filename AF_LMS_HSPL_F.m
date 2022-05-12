% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 2 ******************************************

% This function implements the adaptation of a Hammerstein SAF (HSAF) structure by using the LMS adaptive algorithm
% Where
% F is the adaptive ﬁlter structure;
% x is the input signal sample x[n];
% d is the desired signal sample d[n];
% y is the output signal sample y[n];
% e is the error signal sample e[n].

function [F,y,e] = AF_LMS_HSPL_F(F,x,d)

M = F.M; % Length of the l i n e a r f i l t e r
F.xw(2:M) = F.xw(1:M-1) ; % S h i f t the input delay−l i n e
[F.xw(1),F.af] = ActFunc(x,F.af); % Load a new input i n t o the delay line
for j=2:M % Constructing the matrix U
if F.af.uIndex == F.af.indexes(j)
F.af.gM(j,:) = F.af.gM(j-1,:);
else
F.af.gM(j,:) = zeros(1,4);
end
end
F.af.gM(1,:) = F.af.g; % Load a new vector in matrix U
F.af.indexes(2:M) = F.af.indexes(1:M-1);
F.af.indexes(1) = F.af.uIndex;
y = F.xw.'*F.w; % F i l t e r output
e = d - y; % Error estimation

% LMS weights and c o n t r o l points update −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
F.w = F.w + F.mu*F.xw.*e; % LMS in Eq . ( 1 0 )
if F.af.aftype > 1 % Kind a c t . f . −1 0 1 2 4 5 *
e_av = F.mQ*e;
ii = F.af.uIndex:F.af.uIndex + F.af.P; % P = Spline order
F.af.Q(ii) = F.af.Q(ii) + e_av*( F.w'*F.af.gM).'; % LMS in Eq . ( 1 1 )
end
end
