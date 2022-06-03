
function out=photodetector(out,dt,fo,length,t_total,bitrate,cancel_noise)
    rng('shuffle')
    kbol=1.38e-23; % Boltzmann constant
    TinK=300; % temperature in Kelvin
    Rload=50; % resistance in Ohm
    h=6.626070040e-34; % h Planck
    nQN=cancel_noise*0.5*h*fo/dt; %%%%%%%%%% power spectral density of quantum noise - shot noise
%% Shot noise    
out = out + sqrt(nQN / 2) * (randn(size(out)) + 1i * randn(size(out)));

%% Non linearity and square law
R=0.7 % Responsivity
out=R*abs(out).^2;

%% thermal noise
nTH=cancel_noise*4*kbol*TinK/Rload/(dt);
out=out+sqrt(nTH)*randn(size(out));
%% Bandwidth Limitation
out=butterworth(out, length, t_total, 4,  bitrate, 0); 
out=real(out);  %%DD output

function B=butterworth(A, samples,pulse_duration,N,fc,f0);

dt=pulse_duration/samples;			
j=-(samples/2):1:(samples/2)-1;
f=j*(1/samples);
f=f/dt;
Hpar=1+i*(((f-f0)/fc).^N);
H=Hpar.^(-1);
A=fft(A);
A=fftshift(A);
A=H.*A;
A=ifftshift(A);
A=ifft(A);
B=A;

function A=amplifier(A, gain, fo,dt )
 rng('shuffle')
    nsp=1.5; %%spontaneous emmission
    kbol=1.38e-23; % Boltzmann constant
    TinK=300; % temperature in Kelvin
    Rload=50; % resistance in Ohm
    h=6.626070040e-34; % h Planck
    nQN=1*0.5*h*fo/dt; %%%%%%%%%% power spectral density of quantum noise - shot noise
    
A=A*sqrt(gain);
% Ay=Ay*sqrt(G);
 
n=h*fo*nsp*(gain-1)/dt; %%%%%%%%%% power spectral density of ASE from EDFA

A = A + sqrt(n / 2) * (randn(size(A)) + 1i * randn(size(A)));
% outy = outy + sqrt(n / 2) * (randn(size(outy)) + 1i * randn(size(outy)));
