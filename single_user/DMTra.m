function [gn,en_bar,bn_bar,Nstar,b_bar,SNRdmt]=DMTra(P,NoisePSD,Ex_bar,N,gap)

% Updated J. Cioffi , Spring 2015 - allows complex-channel view
%
% function [gn,en_bar,bn_bar,Nstar,b_bar]=complexdmtra(P,NoisePSD,Ex_bar,N,gap)
%
% CAUTION - when using this program, you should know if your input channel
% P is complex baseband or real in interpretting the outputs, as per
% comments below
%
% INPUTS
% P is the complex baseband pulse response (and v is set equal to length(P)-1)
% Ex_bar is the flat 2-sided power spectral density level so total power divided by
%     positive and negative-frequency-total bandwidth.
% NoisePSD is the flat noise PSD on the same (2-sided) scale as Ex_bar
%    
% N is the DFT size, N>2, and is Nbar from notes for complex channels and N
%     for real channels
% 
%    The sampling rate on each rail (Inphase or Quad) is presumed to be (N+v)/T 
%    where 1/T is DMT symbol rate
% gap is the gap in dB
%
% OUTPUTS
% gn is vector of channel gains
% en_bar is the vector of energy/dim for all subchannels (for complex P
%    double this number to get energy/tone) For real channels, the upper
%    image frequencies will duplicate the lower frequencies, but this is
%    stil energy per real dimension.
% bn_bar is the vector of bit/dim for all subchannels b (for complex P
%    double this number to get bits/tone) For real channels, the upper
%    image frequencies will duplicate lower frequencies, but this is still
%    bits per real dimension for those.
% Nstar is the number of used input-size DFT (real or complex) dimensions 
%     (2*Nstar is number of real dimensions for complex P)
% b_bar is the number of bits per (real) dimension (so 2*(N+v)*b_bar is 
%    total number of bits/symbol for a complex P
% SNRdmt is the equivalet DMT SNR (complex or real)

% dB into normal scale
Noise_var=NoisePSD;
gap=10^(gap/10);
nu=size(P);
nu=max(nu(1),nu(2))-1;

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
% Waterfilling for DMT %
%%%%%%%%%%%%%%%%%%%%%%%

%sort
[gn_sorted, Index]=sort(gn);  % sort gain, and get Index

gn_sorted = fliplr(gn_sorted);% flip left/right to get the largest 
                              % gain in leftside
Index = fliplr(Index);        % also flip index  


num_zero_gn = length(find(gn_sorted == 0)); %number of zero-gain subchannels
Nstar=N - num_zero_gn;    
 	% Number of used channels, 
 	% start from N - (number of zero gain subchannels)

while(1) 
 	K=1/Nstar*(N*Ex_bar+gap*sum(1./gn_sorted(1:Nstar))); 
 	En_min=K-gap/gn_sorted(Nstar);	% En_min occurs in the worst channel
 	if (En_min<0)		
    		Nstar=Nstar-1;  % If negative En, continue with less channels
 	else 
    		break;       % If all En positive, done.
 	end
end

En=K-gap./gn_sorted(1:Nstar); 		% Calculate En
Bn=log2(K*gn_sorted(1:Nstar)/gap); 	% Calculate bn

bn(Index(1:Nstar))=Bn;		% return values in original index
en(Index(1:Nstar))=En;		% return values in original index

		% Since channel is even, need to display 
                                % only half of result
en_bar=en(1:N);
bn_bar=0.5*bn(1:N);

% calculate b_bar
b_bar=1/(N+nu)*(0.5*sum(bn));
SNRdmt=10*log10(gap*(2^(2*b_bar)-1));





