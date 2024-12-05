function error = MSE_fitting_error(All_para,CP_data)

% All_para0 = [EP_real0,EP_imag0,VP_real0,VP_imag0,alpha0];
% 
EP_real = All_para(1,1:12);
EP_imag = All_para(1,13:24);
VP_real = All_para(1,25:36);
VP_imag = All_para(1,37:48);
alpha = All_para(1,49:end);

EP = EP_real + 1i * EP_imag;
VP = VP_real + 1i * VP_imag;

f = 0.05 * [2 3 5 7 11 13 19 23 29 31 37 41];
error = 0;
    for cases = 1:length(alpha)
        prediction = (1 - alpha(cases)) * EP + alpha(cases) * VP;
        error = error + (1/12) * (norm((2 * pi * f) .* (CP_data(cases,:) - prediction),'fro'))^2;
    end
end