function [ax,h]=suplabel(text,whichLabel,supAxes)
% PLaces text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.

% whichLabel is any of 'x', 'y', 'yy', or 't', specifying xlabel, ylabel,
% right-side ylabel, or title

% h = suplabel('text') returns a handle to the text object
% [ax,h]=suplabel('text') returns handles to both the axis and the text
% object

if nargin < 3
    supAxes = [.08 .08 .84 .84];
    ah=findall(gcf,'type','axes');
    if isempty(ah)
        return
    end
    leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
    axBuf=.04;
    set(ah,'units','normalized')
    ah=findall(gcf,'type','axes');
    for k=1:length(ah)
        if strcmp(get(ah(k),'Visible'),'on')
            thisPos=get(ah(k),'Position');
            leftMin=min(leftMin,thisPos(1));
            bottomMin=min(bottomMin,thisPos(2));
            leftMax=max(leftMax,thisPos(1)+thisPos(3));
            bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
        end
    end
    supAxes=[leftMin-axBuf bottomMin-axBuf ...
        leftMax-leftMin+2*axBuf bottomMax-bottomMin+2*axBuf];
end

if ~isstr(text) & (iscell(text) | length(text) > 1)
    for k=1:length(text)
        ax=axes('Position',supAxes);
        set(ax,'visible','off');
        h(k)=text('Position',[.5 1-(k-1)*.1],...
            'Units','normalized',...
            'String',text{k},'Interpreter','none', ...
            'VerticalAlignment','top', ...
            'HorizontalAlignment','center');
    end
    return
end

ax=axes('Position',supAxes);
set(ax,'visible','off');

if strcmp('x',whichLabel)
    set(get(ax,'XLabel'),'Visible','on')
    h=xlabel(text,"FontSize",18);
elseif strcmp('y',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    h=ylabel(text,"FontSize",18);
elseif strcmp('yy',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    h=ylabel(text,"FontSize",18);
    set(h,'Rotation',270)
    set(h,'VerticalAlignment','bottom')
elseif strcmp('t',whichLabel)
    set(get(ax,'Title'),'Visible','on')
    h=title(text,"FontSize",18);
end

if (nargout < 2)
    return
end
