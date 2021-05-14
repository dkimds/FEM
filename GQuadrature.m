function xw = GQuadrature(n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function TriGaussPoints provides the Gaussian points and weights    %
%   for the Gaussian quadrature of order n for the standard triangles.  %
%                                                                       %
%   Input: n   - the order of the Gaussian quadrature (n<=12)           %
%                                                                       %
%   Output: xw - a n by 3 matrix:                                       %
%              1st column gives the x-coordinates of points             %
%              2nd column gives the y-coordinates of points             %
%              3rd column gives the weights                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference
% David Dunavant, High degree efficient symmetrical Gaussian quadrature
% rules for the triangle. International journal for numerical methods in
% engineering. 21(6):1129--1148, 1985.

%--------------------------------------------------------------------------
% Edit)
% The code is modified by switching between 1st column and 2nd column to
% set the reference triangle as (0,0)-(1,0)-(0,1), so the orientation is
% counter-clockwise.
% 
% Edited by Heonkyu Ha
%--------------------------------------------------------------------------

if (n == 3)
    xw=[0.33333333333333    0.33333333333333   -0.56250000000000
        0.20000000000000    0.20000000000000    0.52083333333333
        0.60000000000000    0.20000000000000    0.52083333333333
        0.20000000000000    0.60000000000000    0.52083333333333];
    
elseif (n == 4)
    xw=[0.44594849091597    0.44594849091597    0.22338158967801
        0.10810301816807    0.44594849091597    0.22338158967801
        0.44594849091597    0.10810301816807    0.22338158967801
        0.09157621350977    0.09157621350977    0.10995174365532
        0.81684757298046    0.09157621350977    0.10995174365532
        0.09157621350977    0.81684757298046    0.10995174365532];
    

    
else
    error('Bad input n');
end

end
