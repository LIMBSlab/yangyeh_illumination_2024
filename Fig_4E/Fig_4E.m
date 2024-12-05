for Fish_ID = 4

boostrapped_workspace = load(fullfile(cd, strcat('Fish_0',num2str(Fish_ID),'_new.mat')));

bootstrapped_gain_data = cell(1,size(boostrapped_workspace.beta_Set,2));

bootstrapped_phase_data = cell(1,size(boostrapped_workspace.beta_Set,2));

bootstrapped_gain_fit = cell(1,size(boostrapped_workspace.beta_Set,2));

bootstrapped_phase_fit = cell(1,size(boostrapped_workspace.beta_Set,2));

for cases = 5
figure,

subplot 211
hold on
grid on
for trials = 1:100

bootstrapped_gain_data{1,cases}(trials,:) = abs(boostrapped_workspace.Bootstrapped_CPs(Fish_ID).CP_bootstrapped(cases).CP_bootstrapped_mean{1,trials});
bootstrapped_gain_fit{1,cases}(trials,:) = abs(boostrapped_workspace.CP_model_cell{1,trials}{1,cases});

end

f = [boostrapped_workspace.f, fliplr(boostrapped_workspace.f)];
filled_gain_data = [(mean(bootstrapped_gain_data{1,cases}) + std(bootstrapped_gain_data{1,cases})),fliplr((mean(bootstrapped_gain_data{1,cases}) - std(bootstrapped_gain_data{1,cases})))];
filled_gain_fit = [(mean(bootstrapped_gain_fit{1,cases}) + std(bootstrapped_gain_fit{1,cases})),fliplr((mean(bootstrapped_gain_fit{1,cases}) - std(bootstrapped_gain_fit{1,cases})))];


fill(f,filled_gain_data,[0,0,0]/256,'EdgeColor','none','FaceAlpha',0.2);
fill(f,filled_gain_fit,[255,0,0]/256,'EdgeColor','none','FaceAlpha',0.2);

semilogx(boostrapped_workspace.f, mean(bootstrapped_gain_data{1,cases},'omitnan'),'k','LineWidth',2);
semilogx(boostrapped_workspace.f, mean(bootstrapped_gain_fit{1,cases},'omitnan'),'r','LineWidth',2);

legend('','','data','model')


set(gca,'xScale','log');
set(gca,'yScale','log');
axis([0 2.1 0.01 5]);
ylabel('Bode Gain (cm/cm)')
title('Fig 4E')



subplot 212
hold on
grid on

for trials = 1:100

bootstrapped_phase_data{1,cases}(trials,:) = (rad2deg(unwrap(angle(boostrapped_workspace.Bootstrapped_CPs(Fish_ID).CP_bootstrapped(cases).CP_bootstrapped_mean{1,trials}))));
bootstrapped_phase_fit{1,cases}(trials,:) = (rad2deg(unwrap(angle(boostrapped_workspace.CP_model_cell{1,trials}{1,cases}))));

end

f = [boostrapped_workspace.f, fliplr(boostrapped_workspace.f)];
filled_phase_data = [(mean(bootstrapped_phase_data{1,cases}) + std(bootstrapped_phase_data{1,cases})),fliplr((mean(bootstrapped_phase_data{1,cases}) - std(bootstrapped_phase_data{1,cases})))];
filled_phase_fit = [(mean(bootstrapped_phase_fit{1,cases}) + std(bootstrapped_phase_fit{1,cases})),fliplr((mean(bootstrapped_phase_fit{1,cases}) - std(bootstrapped_phase_fit{1,cases})))];

fill(f,filled_phase_data,[0,0,0]/256,'EdgeColor','none','FaceAlpha',0.2);
fill(f,filled_phase_fit,[255,0,0]/256,'EdgeColor','none','FaceAlpha',0.2);

semilogx(boostrapped_workspace.f, mean(bootstrapped_phase_data{1,cases},'omitnan'),'k','LineWidth',2);
semilogx(boostrapped_workspace.f, mean(bootstrapped_phase_fit{1,cases},'omitnan'),'r','LineWidth',2);

set(gca,'xScale','log');
axis([0 2.1 -250 0]);
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')

end


Luminance = zeros(1,size(boostrapped_workspace.beta_Set,2));
alpha_vector_mat = zeros(100,size(boostrapped_workspace.beta_Set,2));

for idxx = 1:size(boostrapped_workspace.struct,2)
Luminance(idxx) = boostrapped_workspace.struct(idxx).luxMeasured;
end

for index = 1:100
alpha_vector_mat(index,:) = boostrapped_workspace.alpha_vector_cell{index};
end

end
