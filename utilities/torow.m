function x = torow(x)
% x = torow(x)
%
% Converts input vector (of whatever orientation) to be a row vector.
% Faster (built-in) approach: x = x(:)'.
%
% Args:
%   x (numeric, required, positional): a shape-(1, N) or shape-(N, 1)
%     numeric array (vector).
%
% Returns:
%   x (numeric, required, positional): a shape-(1, N) numeric array
%     (row vector).

arguments
  x {mustBeVector(x,'allow-all-empties')}
end

if isempty(x)
  return
end
assert(isvector(x), 'input not a vector')

if size(x, 1) > 1
  x = x';
end