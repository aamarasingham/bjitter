% Matlab code to make Figure 2 for “Spike-centered jitter can mistake temporal structure” (Platkiewicz, Stark, Amarasingham)
% Generates Poisson spike trains and compares "p-value" distributions
% for interval vs. basic (spike-centered) jitter, using a synchrony
% statistic. Synchronies are injected to demonstrate sensitivity differences between basic and interval jitter. 
% Refer to manuscript for details.
% (C) Asohan Amarasingham, 6/5/2016

function []=Appendix_Section_5_3

synch_rate=10; % synchrony process firing in Hz
syn_jitter=.01;  % jitter radius for the synchrony process
			     % synchrony process is generated and then superposed onto the individual spike processes
				 % after injecting a little jitter
frate1=20-synch_rate; % neuron 1 firing in Hz
frate2=20-synch_rate; % neuron 2 firing in Hz

T=1;  % end time in seconds
synch_def=.03;   % spikes x,y synchronous if |x-y|<synch_def in secs
synch_range=[0 1]; % only count synch in this range (ie all synch spikes in neuron 1 \in synch_range)

num_jitter=500;
num_runs=50000;
jitter_width=0.02;
u=rand(num_runs,1);


%%%%%%%

clear pval pvalr pval_int pvalr_int
for ccc=1:num_runs

   orig_syn=0; 
    
    % sample Poisson by sampling exponential ISI's

	% synchrony process
	ISIS_avg=1/synch_rate;  % ISI mean rate
	syn=[ exprnd(ISIS_avg) ];
	while syn(end) < T, syn(end+1)=syn(end) + exprnd(ISIS_avg); end; syn=syn(1:end-1);

    % neuron 1 (base)
    ISI1_avg=1/frate1;  % ISI mean rate
    n1=[ exprnd(ISI1_avg) ];
    while n1(end) < T, n1(end+1)=n1(end) + exprnd(ISI1_avg); end; n1=n1(1:end-1);

    % neuron 2 (base)
    ISI2_avg=1/frate2;  % ISI mean rate
    n2=[ exprnd(ISI2_avg) ];
    while n2(end) < T, n2(end+1)=n2(end) + exprnd(ISI2_avg); end; n2=n2(1:end-1);

	% superpose synchrony process with n1 and n2
	n1=[ n1 syn+syn_jitter*rand(1,length(syn)) ];
	n2=[ n2 syn+syn_jitter*rand(1,length(syn)) ];

    % compute synchrony
    
    orig_syn=synch_compute( n1,n2,synch_def,synch_range );
    orig_synb=orig_syn+.5*rand(1);   % randomized synchrony

    % [basic] jitter, and tabulate synchrony counts

        syn_surr=[]; syn_surrb=[];
        for k=1:num_jitter

            % jitter spikes
            n1_jitt=n1 + (2*jitter_width*(rand(1,length(n1))))-jitter_width;
            n2_jitt=n2;
            
            % compute synchrony
            s=synch_compute( n1_jitt,n2_jitt,synch_def,synch_range );

            syn_surr(k)=s;
            syn_surrb(k)=s+.5*rand(1);   % store synchrony for surrogate j

        end
        
    % [interval] jitter, and tabulate synchrony counts

        syn_surr_int=[]; syn_surrb_int=[];
        for k=1:num_jitter

            % interval jitter (interval length jitter_width*2) spikes for n1
            n1_jitt_int=(jitter_width*2)*floor( n1/(jitter_width*2) ) + (jitter_width*2)*rand( 1,length(n1) );
            %max( n1-n1_jitt )
            n2_jitt=n2;
            
            % compute synchrony
            s=synch_compute( n1_jitt_int,n2_jitt,synch_def,synch_range );

            syn_surr_int(k)=s;
            syn_surrb_int(k)=s+.5*rand(1);   % store synchrony for surrogate j

        end
        
   % compute pvalues
        
   % p(X,R) for basic jitter
   pval(ccc)=(1+sum( syn_surr>=orig_syn))/(num_jitter+1);
   % p_C(X,R) for basic jitter
   pvalr(ccc)=(1+sum( syn_surrb>=orig_synb))/(num_jitter+1);
   % p(X,R) for interval jitter test
   pval_int(ccc)=(1+sum( syn_surr_int>=orig_syn))/(num_jitter+1);
   % p_C(X,R) for randomized interval jitter test
   pvalr_int(ccc)=(1+sum( syn_surrb_int>=orig_synb))/(num_jitter+1);
        


    if mod(ccc,200)==0
        orig_syn,ccc
        binw=.01;
        subplot(2,2,1)
        hold off, cdfplot(pval), title('p(X,R) basic jitter')
        hold on, plot(0:.005:1,0:.005:1,'r-.')
        subplot(2,2,3)
        hold off, cdfplot(pvalr), title('p_C(X,R) basic jitter')  
        hold on, plot(0:.005:1,0:.005:1,'r-.')
        subplot(2,2,2)
        hold off, cdfplot(pval_int), title('p(X,R) interval jitter')
        hold on, plot(0:.005:1,0:.005:1,'r-.') 
        subplot(2,2,4)
        hold off, cdfplot(pvalr_int), title('p_C(X,R) interval jitter') 
        hold on, plot(0:.005:1,0:.005:1,'r-.')
    end
                
end


function synch= synch_compute( n1,n2,synch_def,synch_range );
% computes sychrony between (spike time) vectors n1 and n2  
% sync_def is radius for def of synchrony
% only count synch in the interval specified in synch_range (ie all synch spikes in neuron 1 \in synch_range)
synch=0;
for j=1:length(n1)
    if n1(j)>=synch_range(1) & n1(j)<=synch_range(2) 
        synch=synch+sum( n2>=n1(j)-synch_def & n2<=n1(j)+synch_def );
    end
end            
        
        


