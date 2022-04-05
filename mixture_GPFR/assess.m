function rcc = assess(predc, preds2, testdata, flag, dispdata, savefilename)
%flag = 0: compute rcc only
%flag = 1: compute rcc and plot the results
%flag = 2: compute rcc, plot the results and save the figures


predup = predc + 1.96*sqrt(preds2);    % 95% confidence interval
predlo = predc - 1.96*sqrt(preds2);
yte = testdata{1}(:,end);      %input points 
rcc = zeros(1,2);
rcc(1,1) = rmse(predc, yte);
a = corrcoef(predc, yte);
rcc(1,2) = a(1,2);


if flag > 0    
    tte = testdata{1}(:,1);        %time points 
    
    figure(10)
    clf    
    plot(tte, yte,'b', tte, predc, 'r--', tte, [predup,predlo],'r:', 'LineWidth',1);
    %hold on;
    %plot(dispdata{1}(:,1), dispdata{1}(:,end), 'b', 'LineWidth',1);
    %ylabel('comz');
    
    if flag > 1
        savefilename = ['./figures/' savefilename];
        print('-depsc2', savefilename);
    end
    
    rcc
    disp('Press any key to continue...')
    pause
end