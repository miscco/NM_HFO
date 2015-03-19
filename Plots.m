% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO.cpp C_Column.cpp H_Column.cpp

function Plots(T)

if nargin == 0
    T       	= 120;  		% duration of the simulation
end

mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO.cpp C_Column.cpp H_Column.cpp;

[V_C, V_H] = HFO(T, 0, 0);

L        = max(size(V_C));
timeaxis = linspace(0,T,L);

figure(1)
subplot(211), plot(timeaxis,V_C)
title('Cortex membrane voltage'), xlabel('time in s'), ylabel('V_C in mV')
subplot(212), plot(timeaxis,V_H)
title('HFO membrane voltage'), xlabel('time in s'), ylabel('V_H in mV')

Fs          = L/T;
[Pxx_C,f_C] = pwelch(V_C-mean(V_C), [], [], [], Fs);
n_C         = find(f_C<=200, 1, 'last' );
[Pxx_H,f_H] = pwelch(V_H-mean(V_H), [], [], [], Fs);
n_H         = find(f_H<=200, 1, 'last' );


figure(2)
subplot(211), plot(f_H,log(Pxx_H));
title('HFO power spectrum'), xlabel('Frequency in Hz'), ylabel('Power')
subplot(212), plot(f_C,log(Pxx_C));
title('Cortex power spectrum'), xlabel('Frequency in Hz'), ylabel('Power')