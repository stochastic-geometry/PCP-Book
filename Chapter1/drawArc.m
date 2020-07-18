function drawArc(p, r, th1, th2, varargin)
% number of points
N = ceil(abs((th2 - th1) * 2 / pi * 30));
dTh = (th2 - th1) / N;

X = p(1) + r*cos(th1 + (0:dTh:N*dTh));
Y = p(2) + r*sin(th1 + (0:dTh:N*dTh));
if(~ishold)
    newplot
end
line(X, Y, varargin{:});