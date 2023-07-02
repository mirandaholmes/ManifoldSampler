% plot output from example1, to make sure it's working

% created June 16, 2022

clear;

datafile = 'example1.txt';

xall = load(datafile);
nx = size(xall,1);

ax = 1.5*[-2 3 -2 3];
ifig = 1;

a = zeros(6);
a(1,2) = 1;  % square edges
a(2,3) = 1;
a(3,4) = 1;
a(4,1) = 1;
a(3,5) = 1;  % triangle edges
a(3,6) = 1;
a(5,6) = 1;
a = a+a';

for ix=1:nx
    
    x = xall(ix,:);
    
    plot_framework(x,a,ifig,ax);
    
    if(ix == 1) 
        pause; 
    end
    
    title(['iter = ',num2str(ix)]);
    
end

