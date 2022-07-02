%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x; xPrint = AMfilter(xPhys,'S');
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPrint(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPrint.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPrint.^(penal-1).*ce;
  dv = ones(nely,nelx);
  [xPrint, dc, dv] = AMfilter(xPhys,'S',dc,dv);
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.05;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    xPrint = AMfilter(xPhys,'S');
    if sum(xPrint(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPrint(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPrint); caxis([0 1]); axis equal; axis on; drawnow;
end
end
%


function [ xi, varargout ] = AMfilter( x, baseplate, varargin )
%AMFILTER Applies a virtual additive manufacturing process to a 
%         2D blueprint design input.
%   Possible uses:
%   xi = AMfilter(x)              design transformation, default orientation
%   xi = AMfilter(x, baseplate)   idem, with baseplate orientation specified
%   [xi, df1dx, df2dx,...] = AMfilter(x, baseplate, df1dxi, df2dxi, ...)
%       This includes also the transformation of design sensitivities
% where
%   x : blueprint design (2D array), 0 <= x(i,j) <= 1
%   xi: printed design (2D array)
%   baseplate: character indicating baseplate orientation: 'N','E','S','W'
%              default orientation is 'S'
%              for 'X', the filter is inactive and just returns the input.
%   df1dx, df1dxi etc.:  design sensitivity (2D arrays)

%INTERNAL SETTINGS
P = 40; ep = 1e-4; xi_0 = 0.5; % parameters for smooth max/min functions

%INPUT CHECKS
if nargin==1, baseplate='S'; end 
if baseplate=='X' 
    % bypass option: filter does not modify the blueprint design
    xi = x;
    varargout = varargin;
    return;
end 
nRot=find(upper(baseplate)=='SWNE')-1;
nSens=max(0,nargin-2); 
if nargout~=nSens+1, error('Input/output arguments mismatch.'); end

%ORIENTATION
x=rot90(x,nRot);
xi=zeros(size(x));
for s=1:nSens
    varargin{s}=rot90(varargin{s},nRot);    
end
[nely,nelx]=size(x); 

%AM FILTER =====================
Ns=3;
Q=P+log(Ns)/log(xi_0); 
SHIFT = 100*realmin^(1/P); % small shift to prevent division by 0
BACKSHIFT = 0.95*Ns^(1/Q)*SHIFT^(P/Q);
Xi=zeros(size(x)); keep=zeros(size(x)); sq=zeros(size(x));
% baseline: identity
xi(nely,:)=x(nely,:); % copy base row as-is
for i=(nely-1):-1:1
    % compute maxima of current base row
    cbr = [0, xi(i+1,:), 0] + SHIFT; % pad with zeros
    keep(i,:) = (cbr(1:nelx).^P + cbr(2:(nelx+1)).^P + cbr(3:end).^P);
    Xi(i,:) = keep(i,:).^(1/Q) - BACKSHIFT;
    sq(i,:) = sqrt((x(i,:)-Xi(i,:)).^2 + ep);
    % set row above to supported value using smooth minimum:
    xi(i,:) = 0.5*((x(i,:)+Xi(i,:)) - sq(i,:) + sqrt(ep));
end
%SENSITIVITIES
if nSens
    dfxi=varargin; dfx=varargin; 
    lambda=zeros(nSens,nelx); 
    % from top to base layer:
    for i=1:nely-1
        % smin sensitivity terms
        dsmindx  = .5*(1-(x(i,:)-Xi(i,:))./sq(i,:));
        %dsmindXi = .5*(1+(x(i,:)-Xi(i,:))./sq(i,:)); 
        dsmindXi = 1-dsmindx; 
        % smax sensitivity terms
        cbr = [0, xi(i+1,:), 0] + SHIFT; % pad with zeros
        dmx = zeros(Ns,nelx);
        for j=1:Ns
            dmx(j,:) = (P/Q)*keep(i,:).^(1/Q-1).*cbr((1:nelx)+(j-1)).^(P-1);
        end        
        % rearrange data for quick multiplication:
        qj=repmat([-1 0 1]',nelx,1);
        qi=repmat(1:nelx,3,1); qi=qi(:);
        qj=qj+qi; qs=dmx(:);
        dsmaxdxi=sparse(qi(2:end-1),qj(2:end-1),qs(2:end-1)); 
        for k=1:nSens
            dfx{k}(i,:) = dsmindx.*(dfxi{k}(i,:)+lambda(k,:));
            lambda(k,:)= ((dfxi{k}(i,:)+lambda(k,:)).*dsmindXi)*dsmaxdxi;
        end
    end
    % base layer:
    i=nely;
    for k=1:nSens
        dfx{k}(i,:) = dfxi{k}(i,:)+lambda(k,:);
    end
end

%ORIENTATION
xi=rot90(xi,-nRot);
for s=1:nSens
    varargout{s}=rot90(dfx{s},-nRot);    
end

end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Matthijs Langelaar,                      %
% Department of Precision and Microsystems Engineering,                    %
% Delft University of Technology, Delft, the Netherlands.                  %
% Please sent your comments to: m.langelaar@tudelft.nl                     %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper "An additive manufacturing filter for         %
% topology optimization of print-ready designs", M. Langelaar (2016),      %
% Struct Multidisc Optim, DOI: 10.1007/s00158-016-1522-2.                  %
%                                                                          %
% This code is intended for integration in the 88-line topology            %
% optimization code discussed in the paper                                 %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund,  % 
% Struct Multidisc Optim, 2010, Vol 21, pp. 120--127.                      %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guarantee that the code is   %
% free from errors. Furthermore, the author shall not be liable in any     %
% event caused by the use of the program.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



