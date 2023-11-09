function mustBeLogicalOrListedType(inputVariable, varargin)
% mustBeLogicalOrListedType(inputVariable, variableType)
%
% Function validates if the input variable is logical or a memeber of any
% other specified types.
%
% Args:
%   inputVariable (any, required, positional): A variable of any type (no
%     restrictions).
%   variableType (string or cell of strings, required, positional): a
%     shape-(1, M) character array or a shape-(1, N) cell of character
%     arrays.
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

if ~islogical(inputVariable) && ~ismember(class(inputVariable),cellstr(varargin))
  eid = 'Type:notValid';
  msg = 'Argument must be logical or one of the prescribed types';
  throwAsCaller(MException(eid,msg))
end