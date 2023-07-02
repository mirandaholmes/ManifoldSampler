% plot output from example1, to make sure it's working

% created June 16, 2022

clear;

kpl = 1;

datafile = 'polymer.txt';

xall = load(datafile);
nx = size(xall,1);     % # of data points
n = size(xall,2)/2;    % # of discs


ax = 1.5*sqrt(n)*[-1 2 -1 2];
ifig = 1;

a = zeros(n);
for ii=1:n-1
    a(ii,ii+1) = 1;
    a(ii+1,ii) = 1;
end

for ix=1:kpl:nx
    
    x = xall(ix,:);
    
    plot_framework(x,a,ifig,ax);
    
    if(ix == 1) 
        pause; 
    end
    
    title(['iter = ',num2str(ix)]);
    
end

