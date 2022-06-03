
%% Parameters
 %  fc --> Single Sideband Bandwidth
 %  f0 --> Frequency Detuning from the carrier
 %  Tdelay_in  --> Input Delay (Before the first coupler)
 %  Tdelay2  -->  Feedback Loop Delay 
 % losa & losb  --> Input/ Output coupler  (e.g 0.5 for 50-50)
 % phase   --> Feedback Phase Shifter Value
 
%%
function Output=ROSS(Input, samples,pulse_duration,fc,f0, Tdelay_in, Tdelay2, losa, losb,phase);

dt=pulse_duration/samples;			
j=-(samples/2):1:(samples/2)-1;
f=j*(1/samples);
f=f/dt;
L=1;  % Feedback strength (i.e variable optical attenuator value)
H=0.5*(1+exp(-i*2*pi.*(f-f0)./(fc))); % Mach Zehnder Delayed Interferometer Transfer Function

Input=fft(Input);
Input=fftshift(Input);

Input=H.*Input*sqrt(1-losa.^2).*sqrt(1-losb.^2)./(1+(losa*losb*sqrt(L)).*H.*exp(-i*(2*pi*f*Tdelay2+phase))).*exp(-i*2*pi*f*(Tdelay_in));

Input=ifftshift(Input);
Input=ifft(Input);

Output=Input;

