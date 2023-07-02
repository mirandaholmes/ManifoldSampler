% Verify output from Manifold Sampler, on  trimer

datafile = 'polymer.txt';
pts = load(datafile);
npts = size(pts,1);

nbins = 60;




% get angles
x32 = pts(:,5:6) - pts(:,3:4);
x12 = pts(:,3:4) - pts(:,1:2); %pts(:,1:2) - pts(:,3:4);
cth = sum(x32.*x12,2) ./ sqrt(sum(x32.*x32,2)) ./ sqrt(sum(x12.*x12,2));
th = real(acos(cth));


figure(1)
clf


h = histogram(th,nbins);
edges = h.BinEdges;  % bin edges
c = (edges(1:end-1) + edges(2:end) ) /2; % bin centers
w = edges(2:end) - edges(1:end-1);  % widths

%  theoretical distributions
%tt = ones(size(c));  % flat distribution
k = 7;   % spring constant
tt = exp(k*(1-cos(c)));  % angle spring, in 1-cos(th) ~ th^2/2

% scale theoretical distributions
tt = tt/sum(tt.*w);
tt = tt .*w * npts;

% plot
hold on
plot(c,tt,'*-','Linewidth',2);
hold off


    
    