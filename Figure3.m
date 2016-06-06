% Matlab code to make Figure 3 for “Basic jitter can hallucinate temporal structure” (Platkiewicz, Stark, Amarasingham)
% Illustrates non-uniformity of randomized p-values induced by Poisson approximation.
% (C) Asohan Amarasingham, 6/5/2016

N=500; % Binomial
p=.1;  % parameters
lambda=N*p;  % for Poisson approximation
m=50000;  % number of samples
binw=0.01;
clf

for j=1:m
    
    % Sample Binomial(N,p)
    S = sum( rand(N,1)<=p );

    % compute p-value (poiss approx.)
    pvalp(j)=1-poisscdf(S-1,lambda);
    
    % compute randomized p-value (poiss approx.)
    pvalp_rand(j)=[1-poisscdf(S-1,lambda)] - rand(1)*poisspdf(S,lambda);
    
    % compute p-value (exact binomial)
    pvalb(j)=1-binocdf(S-1,N,p);
    
    % compute randomized p-value (binomial xact)
    pvalb_rand(j)=[1-binocdf(S-1,N,p)] - rand(1)*binopdf(S,N,p);   
    
    %if mod(j,1000)==0, j, end

end

subplot(2,2,1)
histogram(pvalp,0:binw:1,'Normalization','probability')
hold on
plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.') % raw poiss approx. pvalues
title('(a) Raw p-values (Poisson approximation method)')

subplot(2,2,3)
histogram(pvalp_rand,0:binw:1,'Normalization','probability')
hold on
plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.') % randomized poiss approx. pvalues
title('(b) Randomized p-values (Poisson approximation)')


subplot(2,2,2)
histogram(pvalb,0:binw:1,'Normalization','probability')
hold on
plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.') % raw binomial (exact) pvalues
title('(c) Exact p-values (Binomial)')


subplot(2,2,4)
%hist(pvalb_rand,nbins)
histogram(pvalb_rand,0:binw:1,'Normalization','probability')
hold on
plot(0:.005:1,binw*ones( size(0:.005:1)),'r-.') % randomized binomial exact pvalues
title('(d) Randomized exact p-values (Binomial)')
