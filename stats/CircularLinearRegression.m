function  [slope,offset,R_value] = CircularLinearRegression(circ,lin,sign)
% [slope,offset,R_value] = CircularLinearRegression(circ,lin,sign)
%
% Function for fitting regression line to circular-linear data.
%
% Args:
%   circ
%   lin
%   sign
%
% Returns:
%   slope
%   offset
%   R_value
%
% Comments:
%   Modified by Martynas Dervinis (martynas.dervinis@gmail.com)
%
% Authors:
%   Peter Petersen(petersen.peter@gmail.com).

% Parse input
if nargin < 2
  sign = 0;
end

if sign == 0
  a = -0.1:0.0001:0.1;
elseif sign > 0
  a = 0:0.0001:0.1;
elseif sign < 0
  a = -0.1:0.0001:0;
end

% Fit the line
cosinepart=zeros(length(a),1);
sinepart=zeros(length(a),1);
R=zeros(length(a),1);

for i=1:length(a)
  cosinepart(i)=sum(cos(circ(:)-(2*pi*a(i)*lin(:))));
  sinepart(i)=sum(sin(circ(:)-(2*pi*a(i)*lin(:))));
  firstterm=(cosinepart(i)/length(circ))^2;
  secondterm=(sinepart(i)/length(circ))^2;
  R(i)=sqrt(firstterm+secondterm);
end
slope=a(R==max(R));
offset=atan2(sinepart(R==max(R)),cosinepart(R==max(R)));
R_value = max(R);
%figure, plot(a,R);
