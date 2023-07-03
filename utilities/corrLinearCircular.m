function [r, pval] = corrLinearCircular(vec1, vec2, options)
% [r, pval] = corrLinearCircular(vec1, vec2, <type>)
%
% Function performs various correlation analyses on two vectors by demand.
% Vector data can be linear or circular.
%
% Args:
%   vec1 (numeric, required, positional): a shape-(1, N) numeric data
%     array. Values in the data vector must correspond to individual data
%     samples.
%   vec2 (numeric, required, positional): a shape-(1, N) numeric data
%     array. Values in the data vector must correspond to individual data
%     samples.
%   type (char, optional, keyword): a shape-(1, M) character array
%     desribing the type of the correlation analysis to perform on the two
%     data vectors. Can take one of the following values:
%       'Pearson' - compute Pearson's linear correlation parametric
%                   coefficient (type 'help corr' for more info);
%       'Kendall' - compute Kendall's tau (type 'help corr' for more info);
%       'Spearman' - compute Spearman's non-parametric rho (default; type
%                    'help corr' for more info);
%       'circular' - circular correlation coefficient for two circular
%                    random variables (type 'help circ_corrcc' for more
%                    info);
%       'circularnp' - non-parametric circular-circular correlation (type
%                      'help circ_corrccnp' for more info);
%       'circlinear' - correlation coefficient between one circular and one
%                      linear random variables (type 'help circ_corrcl' for
%                      more info);
%       'circlinearnp' - non-parametric circular-linear correlation with
%                        vec1 being circular and vec2 being linear
%                        variables (type 'help circ_corrclnp' for more
%                        info).
%
% Returns:
%   r (numeric): a shape-(1, 1) numeric scalar correlation coefficient.
%   pval (numeric): a shape-(1, 1) numeric scalar p-value associated with
%     the correlation coefficient.
%
% Comments:
%   If you are interested in circular-linear correlation, vec1 is expected
%   to be circular while vec2 is linear.
% 
% Dependencies:
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%   Non-parametric circular statistics functions for Matlab (https://github.com/dervinism/circStatNP).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  vec1 {mustBeVector}
  vec2 {mustBeVector}
  options.type {mustBeMember(options.type,{'Pearson','Kendall','Spearman','circular','circularnp','circlinear','circlinearnp'})} = 'Spearman'
end

% Parse input
vec1 = vec1(:)';
vec2 = vec2(:)';
assert(numel(vec1) == numel(vec2));
idx = ~isnan(vec1) & ~isnan(vec2); % Non-NaN indices
vec1 = vec1(idx);
vec2 = vec2(idx);

% Correlate signals
if strcmpi(options.type,'circular')
  [r, pval] = circ_corrcc(vec1, vec2);
elseif strcmpi(options.type,'circularnp')
  [r, pval] = circ_corrccnp(vec1, vec2);
elseif strcmpi(options.type,'circlinear')
  [r, pval] = circ_corrcl(vec1, vec2);
elseif strcmpi(options.type,'circlinearnp')
  [r, pval] = circ_corrclnp(vec1, vec2);
else
  if sum(vec1) && sum(vec2)
    [r, pval] = corr([vec1' vec2'], 'Type',options.type);
    r = r(2);
    pval = pval(2);
  else
    r = 0;
    pval = 1;
  end
end