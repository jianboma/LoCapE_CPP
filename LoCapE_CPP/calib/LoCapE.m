function [f, eta, XLocap, fL, XLocapall] = LoCapE(x, f1, f2, p, M, L, eta)
% Locap is to perform localized capon algorithm for damped sinusoidal
% signals. This algorithm is to solve the optimization problem
% \min_{h} (Fh)'R(Fh) s.t. (Fh)'Fa(f) = 1; where a(f) is the fourier vector
% at frequency f, and F is the constructed matrix by picking columns from fourier matrix.
%
% Inputs:
%   - x           : input data sequence, a column vector
%   - fc          : frequency component that is interested in
%   - p           : number of columns that will pick at both sides of fc
%   - M           : order of locap filter
%   - L           : points of frequency searching in the range (0, 1]
%   - eta         : a vector containing damping factor
%
% Outputs:
%   - f           : estimated frequency
%   - rho         : estimated damping factor for f
%   - XLocap      : estimated locap energy matrix
%   - fL          : frequency grid of processed range
%   - XLocapall   : estimated locap energy matrix of whole range
%
% use example: [fe, rhoe, ~] = Locap(data, 0.1, 3, 32, 1024);
% References:
%   [1] P. Stoica, R. L. Moses et al., "Spectral analysis of signals," 2004
%   [2] S. Ye, E. Aboutanios, D. S. Thomas, and J. M. Hook, "Localised high resolution
%       spectral estimator for resolving superimposed peaks in nmr signals,"
%       Signal Processing, vol. 130, pp. 343-354, 2017.
%   [3] P. Stoica and T. Sundin, "Nonparametric nmr spectroscopy," Journal of
%       Magnetic Resonance, vol. 152, no. 1, pp. 57-69, 2001.
%
%
% Jianbo Ma <jianboma5@gmail.com>
% EE&T, University of New South Wales

% set default parameter values
N = length(x);
if ( nargin < 2 )
    f1 = -0.5;
end
if ( nargin < 3 )
    f2 = 0.5;
end
if ( nargin < 4 )
    p = 3;
end
if ( nargin < 5 )
    M = 32;
end
if ( nargin < 6 )
    L = 2*N;
end
K = N-M+1;
% K = 20;
if ( nargin < 7 )
    df = 0;
end
% set index vectors
n = 0:N-1;
m = 0:M-1;
k = 0:K-1;
ki = N:-1:N-K+1;%K-1:-1:0;

Tn = 1/L; % frequency interval
pn = floor(L/M); % jump size 
indcol = -p*pn:pn:p*pn; % base index of picked columns

fv = 0:1/L:1-1/L; % frequency grid
% for damping factor
% Le = (1-exp(-2*eta*K))./(1-exp(-2*eta));
% if eta == 0
%     Le = K;
% end


%----------algorithm part-------------%
jL = floor(f1*L)+1;
jH = ceil((f2*L))+1;

%----------------construct signal in a more efficient way
indfreq = jL+ indcol(1):jH + indcol(end); % desired frequencies' index
if indfreq(end)>length(fv)
    indfreq = jL+ indcol(1):length(fv);
end
Lindex = length(indfreq);                 % length of frequencies needed to search


F = exp(-1i*2*pi*fv(indfreq)'*m);

x2 = conj(x(N:-1:1));
XF = x(1:M).';
XB = x2(1:M).';

YF = zeros(Lindex,K); % record signals in frequency domain 
YB = zeros(Lindex,K);
YF(:,1) = F*XF;
YB(:,1) = F*XB;

for kk = 2:K
    YF(:,kk) = (YF(:,kk-1) - x(kk-1)) .* conj(F(:,2))+ x(M+kk-1)*F(:,M); % Maybe related to x2 = [x1[1:end], x[M+1]]
    YB(:,kk) = (YB(:,kk-1) - x2(kk-1)) .* conj(F(:,2)) + x2(M+kk-1)*F(:,M); %this is exactly recursively construct the hankal matrix, that transformed by FM. 
end




XLocap = zeros(length(eta), jH-jL+1);
XLocapall = zeros(length(eta), length(fv));



for neta = 1:length(eta)
    Le = (1-exp(-2*eta(neta)*K))./(1-exp(-2*eta(neta)));
    LeN = (1-exp(-2*eta(neta)*N))./(1-exp(-2*eta(neta))); % for calculating energy
    if eta(neta) == 0
        Le = K;
        LeN = N; % for calculating energy
    end
    for jl = jL:jH
        locap = zeros(length(eta),jH-jL+1);
        
        jindcol = jl+indcol-indfreq(1)+1; % move index of picked columns
        
        a = F(jindcol,:)*exp((-eta(neta)+1i*2*pi*fv(jl))*m.');
        s = (1/Le)*exp((-eta(neta)-1i*2*pi*fv(jl))*(0:K-1)');
        
        X1 = YF(jindcol,:);
        Xi1 = YB(jindcol,:);
        
        R = (1/2)*(X1*X1'+Xi1*Xi1');
        Rinv = inv(R);                                  % inverse of covariance matrix
        
        %   XLocap(neta, jl-jL+1) = N/abs(a'*(Rinv*a))^2*abs(a'*Q*a); % compute energy for given frequency
        XLocap(neta, jl-jL+1) = (a'*Rinv*X1*s/(a'*Rinv*a)); % as in the paper *Le
    end
%     XLocap(:, jl-jL+1) = max(locap(:,jl-jL+1));
%     XLocapall(:,jl-jL+1) = locap(:,jl-jL+1);
end
fL = fv(jL:jH);
if  ( nargout < 1 )
    [YfL, Xrhovec] = meshgrid(fL, rhovec );
    figure;
    mesh(YfL, Xrhovec,10*log10(XLocap));
end

if  ( nargout < 2 )
    val = max(XLocap(:));
    [~, col] = find(ismember(XLocap, val));
    f = fL(col);
end
if  ( nargout < 6 )
    if length(eta) == 1
        [~,locs] = findpeaks(abs(XLocap(1,:)), 'MinPeakHeight',0.05*(max(abs(XLocap(1,:)))));
        f = fL(locs);
        etaE = zeros(1, length(locs));
        XLocapall(:,jL:jH) = XLocap;
    else
        LMaxFinder = vision.LocalMaximaFinder('ThresholdSource','Input port');
        idx = LMaxFinder(abs(XLocap),0.5*max(abs((XLocap(:)))));
        [val, IA, IB] = unique(idx(:,1));
        f = fL(idx(IA,1));
        etaE = eta(idx(IA,2));
        XLocapall(:,jL:jH) = XLocap;
    end
end



