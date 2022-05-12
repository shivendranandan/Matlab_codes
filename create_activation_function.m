% Group 01
% M21ME005- Shivendra Singh
% M21ME009- Shivendra Nandan

%***********************Code 3 ******************************************

% This function creates and initializes nonlinear function implemented by splines
% where:
% AF is the spline nonlinear structure;
% aﬁnit is the type of initialization (linear, Gaussian, random,...). In particular:
% -1: signed sigmoid;
% 0: linear function;
% 1: unsigned sigmoid;
% 2: Gaussian;
% 3: random.
% aftype is the type of spline basis (B, Catmull-Rom, Hermite,...). In particular:
% -1: ﬁxed signed sigmoid;
% 0: ﬁxed linear function;
% 1: ﬁxed unsigned sigmoid;
% 2: ﬁxed Gaussian;
% 3: fexible polynomial function;
% 4: fexible Catmull-Rom spline;
% 5: fexible B-spline;
% 6: ﬂexible Bernstein spline;
% 7: ﬂexible parametric τ-spline;
% 8: ﬂexible Hermite spline;
% 9: ﬂexible Beziér spline;
% 20: ﬂexible quadratic B-spline.
% DeltaX is the space between knots Δx;
% Gain is the x-axis range limits;
% Slope is the initial slope;
% M is the length of the linear ﬁlter M;
% Pord is the polynomial order in case of polynomial nonlinearity Pord

function AF = create_activation_function( afinit, aftype, DeltaX, Gain,Slope, M, Pord)

% Spline a c t i v a t i o n function d e f i n i t i o n and i n i t i a l i z a t i o n −−−−−−−−−−−−−−
if nargin==0 , help create_activation_function; return;end
if nargin < 7
    Pord = 3 ;
    if nargin <6
        M=1;
        if nargin <5
            Slope= 2.16 ;
            if nargin <4
                Gain= 1.1 ;
                if nargin <3
                    DeltaX= 0.4 ;
                    if nargin <2
                        aftype=-1;
                        if nargin <1
                            afinit=-1;
                        end
                    end
                end
            end
        end
    end
end

% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% Check fo n o n l i n e a r i t y type
if afinit == -1
    Slope = Slope*2*(1/Gain);
end

if afinit == 1
    Slope = Slope*4*(1/Gain);
end

if (( afinit<-1) || ( afinit>3))
    fprintf(' Activation function error type\n ');
    aftype = -1;
end

% LUT parameters −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
Table_Length = G1_TabFuncLen(DeltaX,Gain,Slope,afinit);
lut_len = Table_Length; % Length of LUT
s = 0.0; % Linear combiner output
x = 0.0; % Nonlinearity output
uIndex = 1 ; % Span index
uSpline = 0 ; % Local a b s c i s s a

% Polynomial n o n l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 3)
    Q = zeros(Pord,1); % Polynomial n o n l i n e a r i t y
    Q(1) = 1 ;
    g = zeros(1,Pord); % Dot u vector
    gM = zeros(Pord,M); % U matrix ( f o r Hammerstein f i l t e r ) in [ 2 ]
else
    Q = zeros(lut_len,1) ; % Look−Up Table n o n l i n e a r i t y
end
% I f NOT s p l i n e −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype < 4)
    C = 0 ;
end

% Catmul−Rom s p l i n e n o n l l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 4)
    C = 0.5*[-1 3 -3 1 ;
        2 -5 4 -1;
        -1 0 1 0;  % Page 775 , top of column 1
        0 2 0 0];
end

% B s p l i n e n o n l l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 5)
    C = (1/6)*[-1 3 -3 1;
        3 -6 3 0; % Page 774 , bottom of column 2
        -3 0 3 0;
        1 4 1 0 ];
end

% Bernstein polynomial n o n l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 6)
    C = [-1 3 -3 1;
        3 -6 3 0;  % Bernstein polynomials
        -3 3 0 0;
        1 0 0 0];
end

% Parametric s p l i n e n o n l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 7)
    tau = 0.5; % tau = 0 . 5 ==> CR−s p l i n e
    C = [-tau 2-tau tau-2 tau;
        2*tau tau-3 3-2*tau -tau;  % Parametric s p l i n e
        -tau 0 tau 0;
        0 1 0 0];
end

% Hermite s p l i n e n o n l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 8)
    C = [2 -2 1 1;
        -3 3 -2 -1; % Hermite polynomials
        0 0 1 0;
        1 0 0 0];
end

% Bezier s p l i n e n o n l l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 9)
    C = (1/6)*[1 3 -3 1;
        3 -6 3 0;  % Bezier polynomials
        -3 3 0 0;
        1 0 0 0];
end

% Quadratic B−s p l i n e n o n l i n e a r i t y −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if (aftype == 20)
    C = [ 1 1 0;
        -2 2 0;
        1 -2 1];
end
% Spline order −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
P = length(C) - 1 ;

% Dummy arrays f o r l ea r n ing algorithms . See Eqs . ( 5 ) and ( 6 ) −−−−−−−−−−−−
if (aftype > 3)
    g = zeros(1,P+1); % The dot u vector
    gM = zeros(M,P+1); % The U matrix ( f o r Hammerstein SAF) in [ 2 ]
end
indexes = ones(M,1);

% STRUCTURE DEFINITION −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
AF = struct('aftype',aftype,'afinit',afinit,'s',s,'x',x,'Q',Q,'lut_len',lut_len,'uIndex',uIndex,'uSpline',uSpline,'Slope',Slope,'Gain',Gain,'DeltaX',DeltaX,'P',P,'Pord',Pord,'C',C,'g',g,'gM',gM,'indexes',indexes);
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

% For s p l i n e i n t e r p o l a t i o n −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
if ( aftype>1 && aftype ~= 3)
    LutSlope = (Table_Length - 1)/2.0; % New slope
    X = -LutSlope*DeltaX;
    for j=1 : Table_Length % Table_Length
        AF.Q( j) = G1_FUNC( X, Gain, Slope, afinit) ;
        X = X + DeltaX;
    end
end
end
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

