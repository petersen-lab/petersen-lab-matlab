function [radius, rotationCentre] = effectiveRotationRadius(pcaScore, options)
% [radius, rotationCentre] = effectiveRotationRadius(pcaScore, <nPCs>)
%
% Function calculates the effective rotation radius for a single roation in
%   the principal component space. Effective here means the mean of all
%   distances between all PC score sample points and the centre of the
%   rotation.
%
% Args:
%   pcaScore (numeric, required, positional): a shape-(N, M) numeric array
%     of principal component analysis output scores with columns
%     corresponding to individual principal components and rows
%     corresponding to observations. The number of rows should be limited
%     approximately to a single rotation unless one is interested in
%     periods longer than a single rotation.
%   nPCs (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the number of components to use for calculating the radius
%     (default=2).
%
% Returns:
%   radius (numeric): a shape-(1, 1) numeric scalar corresponding to the
%     rotation radius.
%   rotationCentre (numeric): a shape-(1, <nPC>) numeric array
%     corresponding to the rotation centre in the <nPC>-dimensional PC
%     space.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  pcaScore (:,:) {mustBeNumeric}
  options.nPCs (1,1) {mustBeNumeric,mustBePositive} = 2
end

rotationCentre = repmat(mean(pcaScore(:,1:options.nPCs)), ...
  size(pcaScore(:,1:options.nPCs),1),1);
rotationRadii = vecnorm((pcaScore(:,1:options.nPCs) - rotationCentre)');
radius = mean(rotationRadii);
rotationCentre = rotationCentre(1,:);