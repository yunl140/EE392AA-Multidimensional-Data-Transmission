function [gn,en,bn,Nstar,b_bar_check,margin]=DMTma(P,NoisePSD,Ex_bar,b_bar,N,gap)

% EE379C 2015 Spring,
% J. Cioffi - allows complex-channel view
%
% CAUTION - when using this program, you should know if your input channel
% is P is complex baseband or real in interpretting the outputs as per
% comments below

% INPUTS
% P - the complex or real baseband pulse response; user must know which.
%     v=length(P) -1
% NoisePSD and Ex_bar - noise power and energy/dimension or cpx-sample. 
%    These two parameters need to be consistent in terms of # dimensions.
% b_bar is the target normalized bit rate in bits/real dimension (so 1/2
%     sum of bits/tone for complex P)
% N is the DFT size, N>2, and is Nbar from notes for complex channels and N
%     for real channels
% The sampling rate (on each of Inphase or Quad for complex) is (N+v)/T 
%    where 1/T is DMT symbol rate
% gap is the gap in dB
%
% OUTPUTS
% gn is vector of channel gains or magnitude-squared |P|^2 / NoisePSD
% en - vector of energy allocation. 
%      Per real dimension for real channels; 
%      Per tone for complex channel.
%    For real channels, the upper image frequencies duplicate the lower 
%      frequencies.
%    For complex channels, each En output is per that tone (= cmplx dim)
% bn is the vector of bit/dim for all subchannels 
%    For complex channels, double bn_bar for bits/tone. 
%    For real channels, the upper image frequencies duplicate lower 
%      frequencies, but bn_bar remains bits per real dimension.
% Nstar is the number of used input-size DFT (real or complex) dimensions 
%     (2*Nstar is number of real dimensions for complex chan)
% b_bar is the number of bits per (real) dimension 
%     (so 2*(N+v)*b_bar is total number of bits/symbol for a complex chan )
% SNRdmt is the equivalent DMT "geometric" SNR (complex or real) in dB

Noise_var=NoisePSD;
gap=10^(gap/10);
nu=length(P)-1;

% initialization
en=zeros(1,N);
bn=zeros(1,N);
gn=zeros(1,N);
Hn = zeros(1,N);

% subchannel center frequencies
f=0:1/N:1-1/N;

% find Hn vector
for i=1:length(P)
	Hn=Hn+P(i)*exp(j*2*pi*f*(i-1)); 
end

% find gn vector
gn=abs(Hn).^2/Noise_var;

%%%%%%%%%%%%%%%%%%%%%%%
% MA waterfilling %
%%%%%%%%%%%%%%%%%%%%%%%


%sort
[gn_sorted, Index]=sort(gn);  % sort gain, and get Index

gn_sorted = fliplr(gn_sorted);% flip left/right to get the largest 
                              % gain in leftside
Index = fliplr(Index);        % also flip index  

num_zero_gn = length(find(gn_sorted == 0)); %number of zero gain subchannels
Nstar=N - num_zero_gn;    
 	% Number of used channels, 
 	% start from N - (number of zero gain subchannels)
b=(N+nu)*b_bar;
Klog2tilde=2*b-(1/log(2))*sum(log(gn_sorted(1:Nstar)));

while(1) 
 	Klog2=(1/Nstar)*Klog2tilde + log(gap)/log(2);
 	En_min=2^(Klog2)-gap/gn_sorted(Nstar);	% En_min occurs in the worst channel
 	if (En_min<0)
            Klog2tilde=Klog2tilde+(1/log(2))*log(gn_sorted(Nstar));
    		Nstar=Nstar-1;  % If negative En, continue with less channels
 	else 
    		break;       % If all En positive, done.
 	end
end
K=2^(Klog2);
En=K-gap./gn_sorted(1:Nstar); 		% Calculate En
Bn=.5*log2(K*gn_sorted(1:Nstar)/gap); 	% Calculate bn

bn(Index(1:Nstar))=Bn;		% return values in original index
en(Index(1:Nstar))=En;		% return values in original index

en_bar=en;
bn_bar=bn;

% calculate b_bar
b_bar_check=1/(N+nu)*(sum(bn));

%calculate margin
margin=10*log10(N*Ex_bar/sum(en));


