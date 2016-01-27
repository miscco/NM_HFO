% mex command is given by: 
% mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO_mex.cpp CA3_Column.cpp Cortical_Column.cpp

function Plots(T)

if nargin == 0
    T       	= 120;  		% duration of the simulation
end

mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O3" HFO_mex.cpp CA3_Column.cpp Cortical_Column.cpp;

[V_C, V_H, Y_H] = HFO_mex(T, 0, 0);

L        = max(size(V_C));
timeaxis = linspace(0,T,L);

figure(1)
subplot(311), plot(timeaxis,V_C)
title('Cortex membrane voltage'), xlabel('time in s'), ylabel('V_C in mV')
subplot(312), plot(timeaxis,V_H)
title('HFO membrane voltage'), xlabel('time in s'), ylabel('V_H in mV')
subplot(313), plot(timeaxis,V_H)
title('HFO synaptic drive'), xlabel('time in s'), ylabel('Y_H in mV')

Fs          = L/T;
[Pxx_C,f_C] = pwelch(V_C-mean(V_C), [], [], [], Fs);
n_C         = find(f_C<=400, 1, 'last' );
[Pxx_H,f_H] = pwelch(V_H-mean(V_H), [], [], [], Fs);
n_H         = find(f_H<=400, 1, 'last' );
[Pxx_Y,f_Y] = pwelch(Y_H-mean(Y_H), [], [], [], Fs);
n_Y         = find(f_Y<=400, 1, 'last' );

figure(2)
subplot(311), plot(f_C(1:n_C),log(Pxx_C(1:n_C)));
title('Cortex power spectrum'), xlabel('Frequency in Hz'), ylabel('Power')
subplot(312), plot(f_H(1:n_H),log(Pxx_H(1:n_H)));
title('HFO power spectrum'), xlabel('Frequency in Hz'), ylabel('Power')
subplot(313), plot(f_Y(1:n_Y),log(Pxx_Y(1:n_Y)));
title('HFO drive power spectrum'), xlabel('Frequency in Hz'), ylabel('Power')