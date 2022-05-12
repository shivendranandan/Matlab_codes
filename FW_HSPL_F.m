% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 5 ******************************************

% This function evaluates the output of a Hammerstein nonlinear structure implemented by a HSAF architecture
% where
%  F is the adaptive ﬁlter structure;
% x is the input signal sample x[n];
% y is the output signal sample y[n];
% s is the linear combiner input array sn.
function [F,y,s] = FW_HSPL_F(F,x)

M = F.M; % Length of the l i n e a r f i l t e r

F.xw(2:M) = F.xw(1:M-1) ; % S h i f t of the input delay−l i n e
F.xw(1) = x; % Load a new input i n t o the delay−l i n e

s = zeros (M,1); % Linear combiner input array
for i=1:M
[s(i),F.af] = ActFunc(F.xw(i),F.af); % Evaluating the n o n l i n e a r i t y output
end
y = s'.*F.w; % F i l t e r output
end