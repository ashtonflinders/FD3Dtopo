function cMap=cmap_div(clim,cmax)
%
% cmap_div: Generate colormap
%
% Usage: cMap=cmap_div(clim,cmax)
%

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Added help information, but uncomplete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default parameters
if nargin < 2, clim=0;cmax=0; end

%ct=hex2rgb('BF0000');
ct=hex2rgb('660000');
%ct=hex2rgb('4D0000');

%cb=hex2rgb('2B0080');
%cb=hex2rgb('101073');
%cb=hex2rgb('260073');
%cb=hex2rgb('300073');
%cb=hex2rgb('001D73');
cb=hex2rgb('001A66');

C=jetwr; C=flipud(C); njet=size(C,1);

cMap=C; ct0=C(end,:); cb0=C(1,:);
if clim>0
   ncmp=ceil( (cmax/clim-1)*njet/2 );
   cMap=[makeColorMap(cb,cb0,ncmp); cMap];
   cMap=[cMap; makeColorMap(ct0,ct,ncmp)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       jewtr colormap                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = jetwr(m,m0);
%JET	Variant of HSV.
%	JET(M), a variant of HSV(M), is the colormap used with the
%	NCSA fluid jet image.
%	JET, by itself, is the same length as the current colormap.
%	Use COLORMAP(JET).
%
%	See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%	modified from jet by Rafi Katzman
%
% gir is the amount of green in the red side:  0 <= gir <= 1
% rd  is the darkenss of the red:   0 <= rd <= 1
gir = 0.75;
rd = 0.5;
%
if nargin < 2, m0=0; end
if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = max(round((m-m0)/5),1);
x = (1:n)'/n;
x2 = gir * x;
x3 = gir + (1 - gir) * x;
y = (n/2:n)'/n;
y2 = 1 - rd * (1 - y);
z = (n/2:n)'/n;
e1 = ones(length(x),1);
e2 = ones(length(z),1);
e3 = ones(m0,1);
r = [y2; e1; e1; e3; flipud(x); 0*e1; 0*y];
g = [0*y; x2; x3; e1; e3; flipud(x); 0*y];
b = [0*y; 0*e1; x; e1; e1; e3; flipud(y)];
J = [r g b];
size(J);
while size(J,1) > m
   J(1,:) = [];
   if size(J,1) > m, J(size(J,1),:) = []; end
end
