function [equationCoeffs, parameters, RSquare, ellipse] = fitPlotEllipse(xy, options)
% [parameters, RSquare] = fitPlotEllipse(xy, type)
%
% Function fits an ellipse to data and displays the fit.
%
% Args:
%   xy (numeric, required, positional): a shape-(N, 2) numeric array of
%     data points in x and y coordinates. The first column of the matrix
%     corresponds to the x coordinates while the second one corresponds to
%     the y coordinates.
%   type (char, optional, keyword): a shape-(1, M) character array
%     describing fit type. Two fit types are available:
%     'ellipseHyperbola' - uses EllipseFitByTaubin function to fit the data
%                          (default). It returns an ellipse or a hyperbola,
%                          if the latter approximates the data better.
%     'ellipseOnly' - uses EllipseDirectFit to fit the data and always
%                     returns an ellipse.
%   draw (logical, optional, keyword): a shape-(1, 1) logical scalar
%     indicating whether to draw the fitted ellipse (default=false).
%
% Returns:
%   equationCoeffs (numeric): a shape-(6, 1) numeric array of algebraic
%     coefficients of the fitted ellipse equation in the form
%     [a; b; c; d; e; f], so that ax^2 + bxy + cy^2 +dx + ey + f = 0.
%   parameters (struct): a shape-(1, 1) structure scalar with the geometric
%     ellipse parameters with fields radii, U, and x0, where y = U'*(x-x0),
%     z = y./radii, and |z|^2 = sum(z.^2) = 1.
%   RSquare (numeric): a shape-(1, 1) numeric scalar providing an estimate
%     of explained variance by the fit (R-square).
%   ellipse (numeric): a shape-(2, 181) numeric array of x and y
%     coordinates of points on the fitted ellipse that if plotted would
%     show the fitted ellipse.
%
% Comments:
%   This function uses EllipsePrj/EllAlg2Geo function to convert from the
%   algebraic to the polar representations.
%
% Dependencies:
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   EllipseFitByTaubin
%     (https://se.mathworks.com/matlabcentral/fileexchange/22683-ellipse-fit-taubin-method).
%   EllipseDirectFit
%     (https://se.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method).
%   EllipsePrj
%     (https://se.mathworks.com/matlabcentral/fileexchange/27711-euclidian-projection-on-ellipsoid-and-conic).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com

arguments
  xy (:,2) {mustBeNumeric}
  options.type (1,:) {mustBeMember(options.type,{'ellipseHyperbola','ellipseOnly'})} = 'ellipseHyperbola'
  options.draw (1,1) {logical} = false
end

% Fit the ellipse to data
if strcmpi(options.type, 'ellipseHyperbola')
  equationCoeffs = EllipseFitByTaubin(xy);
elseif strcmpi(options.type, 'ellipseOnly')
  equationCoeffs = EllipseDirectFit(xy);
end

% Convert the equation coefficients to the polar (parametric) repsentation
A = equationCoeffs(1);
B = equationCoeffs(2);
C = equationCoeffs(3);
D = equationCoeffs(4);
E = equationCoeffs(5);
F = equationCoeffs(6);
H = [A,   B/2;
     B/2, D];
g = 0.5*[C; E];
c = F;
[parameters.radii, parameters.U, parameters.x0] = EllAlg2Geo(H, g, c);

% Points on ellipse, parametric 
Ea = parameters.U*diag(parameters.radii);
theta = linspace(0,2*pi,181);
ellipse = parameters.x0 + Ea*[cos(theta); sin(theta)];

% Implicit function on grid
%[minxy, maxxy] = bounds(ellipse,2);
%x = linspace(minxy(1),maxxy(1));
%y = linspace(minxy(2),maxxy(2));
%[X,Y] = meshgrid(x,y);
%XY = [X(:) Y(:)]';
%Z = reshape(sum(XY.*(H*XY + g),1) + c, size(X)); % == (A*x^2)+(B*x*y)+(C*x)+(D*y^2)+(E*y)+F
%Z = reshape(Z, size(X));

% Draw the ellipse
if options.draw
  figure
  %imagesc(x,y,Z);
  %hold on
  %contour(x,y,Z,[0 0],'ro');
  plot(ellipse(1,:),ellipse(2,:),'k');
  axis equal;
end

% Goodness of fit estimation
residualSumSquare = 0;
for iResidual = 1:size(xy,1)
  residuals = vecnorm([ellipse(1,:)' - xy(iResidual,1) ellipse(2,:)' - xy(iResidual,2)]');
  residualSumSquare = residualSumSquare + min(residuals)^2;
end
errorSumSquare = sum(vecnorm([xy(:,1) - mean(xy(:,1)) xy(:,2) - mean(xy(:,2))]').^2);
RSquare = max([0 1 - (residualSumSquare/errorSumSquare)]);