% plot framework and given flex, stress, in 2d

% created jan 26 2018


function plot_framework(points,a,whichfig,ax)

n = size(a,1);


% color map for stress 0.2 for final plot, 0.7 for now
lcolr = 0.7*[1 1 1];    % line color
lstyl = '-';
lw = 3;
msize = 8;   % marker size for points
iftext = 1;
if(n > 200)
    lw = 2;
    msize = 6;
end
if(n > 400)
    msize = 2;
    iftext = 0;
    lw = 1;
end


cmap = lines(7);  % default color map
cblue1 = [0    0.4470    0.7410];
cblue2 = [0.6    0.9    1];
cblue3 = [0.2    0.9    0.7];
corange = [0.8500    0.3250    0.0980];
cred = [0.6350    0.0780    0.1840];
cpurple = [0.5,0,0.5];


% -----   Construct needed variables   -----

x = points(1:2:end);
y = points(2:2:end);
[rr,cc] = find(triu(a));
edges = [rr,cc];
nb = length(rr);  % number of regular edges


% -----  Set up plot  -----

figure(whichfig);
clf
hold on

% -----   Plot edges and stress  -----

% Plot Regular edges
i1 = edges(1:nb,1);
i2 = edges(1:nb,2);
x1 = x(i1)'; x2 = x(i2)';
y1 = y(i1)'; y2 = y(i2)';
ih = plot([x1,x2]',[y1,y2]','Linewidth',lw,'LineStyle',lstyl,'Color',lcolr);


% -----   Plot points   -----
% Plot interior points
plot(x,y,'o','Markersize',msize,'MarkerEdgeColor','k','MarkerFaceColor',cblue1);



% -----   Plot text   -----
if(iftext)
    tshift = 0.04;
    for jx=1:length(x)
        text(x(jx)+tshift,y(jx),num2str(jx),'HorizontalAlignment','left')
    end
end


%axis equal

axis(ax);
axis square
hold off






