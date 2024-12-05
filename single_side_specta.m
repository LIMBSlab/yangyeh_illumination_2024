function [P1,f] = single_side_specta(X,mean_subtract,Fs,fig,xlim_vec)


X(isnan(X))=[];
if mean_subtract
    X = X - mean(X);
end

Y = fft(X);
L = length(X);

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.

P2 = abs(Y/L);
P1 = P2(1:(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
% Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

f = (Fs/L)*(0:(L/2));

if fig

    figure('color','white')
    set(gca,'linewidth',1,'fontsize',14)
    hold on
    stem(f,P1,"filled",'LineWidth',1.5) 
    title('Single-Sided Amplitude Spectrum of X(t)','FontWeight','n')
    xlabel('f (Hz)','FontSize',14)
    ylabel('|P1(f)|','FontSize',14)
    xlim(xlim_vec)
end