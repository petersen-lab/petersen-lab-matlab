% This program bring all the figures to the front
% Liming Shi 2016-11-20 in Aalborg, Denmark
% How to use it? --> click-Home-set Path-Add Folder (crate/use a customized folder
% for all of your customized functions)-put this program inside it.
% The ploting sequence is from left to right and low to high. After the
% screen is filled, ploting starts from high to low but still left to right.
function bring_figure_front()
handle_set = findobj('Type', 'Figure');
handle_set_num=[];
for cnt=1:length(handle_set)
    handle_set_num=[handle_set_num handle_set(cnt).Number];
end
handle_set_num=sort(handle_set_num);

x_len=560/(6/5);
y_len=420/(6/5);
x_position=0;
y_position=0;
Pix_SS = get(0,'screensize');
for count=1:length(handle_set_num)
    H1=sfigure(handle_set_num(count));
    set(H1,'Position',[x_position,y_position,x_len,y_len]);
    format_plot_ls;
    if x_position+x_len<Pix_SS(3) 
        x_position=x_position+(x_len-1/15*x_len);
    else
        x_position=0;
        if y_position==0 || y_position+2*y_len<Pix_SS(4)
            y_position=y_position+y_len;
        else
            y_position=y_position-y_len/2;
        end
    end
end
for count=1:length(handle_set_num)
    figure(handle_set_num(count));
end
return
function format_plot_ls()
    set(gca, 'Box', 'off')
    % set(gca, 'PlotBoxAspectRatio', [2 1 1])
    set(gca, 'TickDir', 'out')
    set(gca, 'XColor', [0.5 0.5 0.5]);
    set(gca, 'YColor', [0.5 0.5 0.5]);
    % set(gca, 'Color', [255/255 255/255 220/255]); % light yellow
    set(gca, 'Color', [255/255 255/255 255/255]); % white
    set(gca, 'LineWidth', 0.01)
    set(gca, 'Visible', 'on')
    set(gca,'FontSize',12)
return
function h = sfigure(h)
% In Windows, MATLAB R14's FIGURE behaviour can be annoying.
% Namely, whenever we invoke figure/figure(X), it immediately jumps to the foreground, 
% and steals focus from whatever other window was active,
% preventing you from switching to other processes, the MATLAB editor window, etc... 
% (if you are repeatedly switching between figures/displaying results in a loop).
% This wrapper to figure makes FIGURE silent -- ie. if you had figure X in the background,
% then the call figure(X) will keep that figure in the background.

% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

if nargin>=1 
	if ishandle(h)
		set(0, 'CurrentFigure', h);
	else
		h = figure(h);
	end
else
	h = figure;
end
return

