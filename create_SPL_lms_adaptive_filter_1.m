% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 4 ******************************************


% This function creates and initializes the structure implementing a WSAF architecture and adapted by the LMS adaptive algorithm.
% Where
% lms_af is the adaptive ﬁlter structure;
% M is the length of the linear ﬁlter M;
% mu is the step-size of the linear ﬁlter µ;
% mQ is the step-size of the nonlinearity µQ;
% delta is the regularization parameter δ;
% af is the spline nonlinearity structure.


function lms_af = create_SPL_lms_adaptive_filter_1( M, mu, mQ, delta, af)

% LMS SAF adaptive f i l t e r d e f i n i t i o n and i n i t i a l i z a t i o n −−−−−−−−−−−−−−−−−
w = zeros(M,1); % Linear f i l t e r taps
w(floor(M/2)) = 1 ; % Adaptive f i l t e r weights i . c .
dI = delta; % Regularizing parameter
xd = zeros(M,1) ; % Buffer of the desired s i g n a l
xw = zeros(M,1) ; % Buffer of the f i l t e r s t a t u s
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
lms_af = struct('M',M,'mu',mu,'mQ',mQ,'dI',dI,'w',w,'xd',xd,'xw',xw,'af',af);
end
