% код по загрузке промежуточных данных, сохранённых по способу ОВ
% 17.12.2020

switch ind
    case 1
        data3=load('From Lisbon STS\Output\benchmark_O2_Vib_noWall_1Torr_O2X_1\ChemDensitiesTime.txt');
    case 2
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_1Torr_O2X_0p7\ChemDensitiesTime.txt');
    case 3
        data3=load('From Lisbon STS\Output\benchmark_O2_Vib_noWall_10Torr_O2X_1\ChemDensitiesTime.txt');
    case 4
        data3=load('From Lisbon\Output\benchmark_O2_noVib_noWall_10Torr_O2X_0p7\ChemDensitiesTime.txt');
    otherwise
        data3=zeros(1,12);
end

use_Y=1;
if use_Y
    win_1=2; win_2=3;
    figure('Position', [70 0 1200 610])
    subplot(win_1, win_2, 1)
    dddraw(data3(:,1), data3(:,5), X, Y(:,num_O2));
    title('nO2(X, 0)')
    subplot(win_1, win_2, 2)
    dddraw(data3(:,1), data3(:,4), X, ...
                    sum(Y(:,num_O2:num_O2-1+O2.num_vibr_levels(1)), 2));
    title('nO2(X, all)')
    subplot(win_1, win_2, 3)
    dddraw(data3(:,1), data3(:,2), X, Y(:,num_O2a));
    title('nO2(a)')
    subplot(win_1, win_2, 4)
    dddraw(data3(:,1), data3(:,3), X, Y(:,num_O2b));
    title('nO2(b)')
    subplot(win_1, win_2, 5)
    dddraw(data3(:,1), data3(:,end-5), X, Y(:,num_O));
    title('nO(3P)')
    subplot(win_1, win_2, 6)
    dddraw(data3(:,1), data3(:,end-4), X, Y(:,num_O1D));
    title('nO(1D)')

    figure('Position', [70 0 1200 610])
    win_1=2; win_2=3;
    subplot(win_1, win_2, 1)
    dddraw(data3(:,1), data3(:,end-1), X, Y(:,num_O3));
    title('nO3')
    subplot(win_1, win_2, 2)
    dddraw(data3(:,1), data3(:,end-6), X, Y(:,num_O2p));
    title('nO2+')
    subplot(win_1, win_2, 3)
    dddraw(data3(:,1), data3(:,end-3), X, Y(:,num_Op));
    title('nO+')
    subplot(win_1, win_2, 4)
    dddraw(data3(:,1), data3(:,end-2), X, Y(:,num_Om));
    title('nO-')
    subplot(win_1, win_2, 5)
    dddraw(data3(:,1), data3(:,end), X, Y(:,num_O3exc));
    title('nO3(exc)')
    %%
    subplot(win_1, win_2, 6)
    % loglog(tm(tm_i), print_data_test(tm_i,1+1), 'LineWidth', 1.2, ...
    %                                                 'color', [1 0.7 0])
    % hold on
    % for i=2:O2.num_vibr_levels(1)
    %     loglog(tm(tm_i), print_data_test(tm_i,1+i), 'LineWidth', 1.2, ...
    %     'color', [1 0.7*(1-i/O2.num_vibr_levels(1)) i/O2.num_vibr_levels(1)])
    % end
    title('nO2(X, i)')
    loglog(X, Y(:,1), 'LineWidth', 1.2, ...
                                                    'color', [1 0.7 0])
    hold on
    for i=2:O2.num_vibr_levels(1)
        loglog(X, Y(:,i), 'LineWidth', 1.2, ...
        'color', [1 0.7*(1-i/O2.num_vibr_levels(1)) i/O2.num_vibr_levels(1)])
    end
%     ylim([1e-10 1e25])
    %%
    % figure('Position', [70 0 1200 610])
    % semilogy(0:41, print_data_test(end,2:1+O2.num_vibr_levels(1)))
    if min(Y(:,1:end), [], 'all')<0
        warning('At least one value of n_ci lower than 0')
    end
else
load print_data_test.txt
% print_data_test=print_data_test1;
tm = print_data_test(:,1); % MP - массив загруженных данных, tm - вектор с шагами интегрирования
    tm_i = [1]; % массив с значащими индексами
    for i = 2:length(tm)-1
        ti = tm(i); tip = tm(i+1);
        if (ti == tip) || (ti==0)
            continue;
        else
            tm_i = [tm_i; i];
            while tm(tm_i(end))>tip
                tm_i=tm_i(1:end-1);
            end
        end
    end

win_1=2; win_2=3;
figure('Position', [70 0 1200 610])
subplot(win_1, win_2, 1)
dddraw(data3(:,1), data3(:,5), tm(tm_i), print_data_test(tm_i,1+num_O2));
title('nO2(X, 0)')
subplot(win_1, win_2, 2)
dddraw(data3(:,1), data3(:,4), tm(tm_i), ...
 sum(print_data_test(tm_i,1+num_O2:1+num_O2-1+O2.num_vibr_levels(1)), 2));
title('nO2(X, all)')
subplot(win_1, win_2, 3)
dddraw(data3(:,1), data3(:,2), tm(tm_i), print_data_test(tm_i,1+num_O2a));
title('nO2(a)')
subplot(win_1, win_2, 4)
dddraw(data3(:,1), data3(:,3), tm(tm_i), print_data_test(tm_i,1+num_O2b));
title('nO2(b)')
subplot(win_1, win_2, 5)
dddraw(data3(:,1), data3(:,end-5), tm(tm_i), print_data_test(tm_i,1+num_O));
title('nO(3P)')
subplot(win_1, win_2, 6)
dddraw(data3(:,1), data3(:,end-4), tm(tm_i), print_data_test(tm_i,1+num_O1D));
title('nO(1D)')

figure('Position', [70 0 1200 610])
win_1=2; win_2=3;
subplot(win_1, win_2, 1)
dddraw(data3(:,1), data3(:,end-1), tm(tm_i), print_data_test(tm_i,1+num_O3));
title('nO3')
subplot(win_1, win_2, 2)
dddraw(data3(:,1), data3(:,end-6), tm(tm_i), print_data_test(tm_i,1+num_O2p));
title('nO2+')
subplot(win_1, win_2, 3)
dddraw(data3(:,1), data3(:,end-3), tm(tm_i), print_data_test(tm_i,1+num_Op));
title('nO+')
subplot(win_1, win_2, 4)
dddraw(data3(:,1), data3(:,end-2), tm(tm_i), print_data_test(tm_i,1+num_Om));
title('nO-')
subplot(win_1, win_2, 5)
dddraw(data3(:,1), data3(:,end), tm(tm_i), print_data_test(tm_i,1+num_O3exc));
title('nO3(exc)')
%%
subplot(win_1, win_2, 6)
% loglog(tm(tm_i), print_data_test(tm_i,1+1), 'LineWidth', 1.2, ...
%                                                 'color', [1 0.7 0])
% hold on
% for i=2:O2.num_vibr_levels(1)
%     loglog(tm(tm_i), print_data_test(tm_i,1+i), 'LineWidth', 1.2, ...
%     'color', [1 0.7*(1-i/O2.num_vibr_levels(1)) i/O2.num_vibr_levels(1)])
% end
title('nO2(X, i)')
loglog(X, Y(:,1), 'LineWidth', 1.2, ...
                                                'color', [1 0.7 0])
hold on
for i=2:O2.num_vibr_levels(1)
    loglog(X, Y(:,i), 'LineWidth', 1.2, ...
    'color', [1 0.7*(1-i/O2.num_vibr_levels(1)) i/O2.num_vibr_levels(1)])
end
ylim([1e-10 1e25])
%%
% figure('Position', [70 0 1200 610])
% semilogy(0:41, print_data_test(end,2:1+O2.num_vibr_levels(1)))
if min(print_data_test(tm_i(:),2:end), [], 'all')<0
    error('At least one value of n_ci lower than 0')
end
end
%%
if ind<5
    figure
%     for i=0:41
% %     dddraw(data3(:,1), data3(:,2), tm(tm_i), print_data_test(tm_i,1+num_O2));
%     semilogx(data3(:,1), data3(:,5+i), ...'LineWidth', 1.5,...
%         tm(tm_i), print_data_test(tm_i,1+num_O2+i), 'LineWidth', 1.2)
% %     pause
%     end
loglog(data3(:,1), data3(:,5), 'LineWidth', 1.2, ...
                                                'color', [1 0.7 0])
hold on
loglog(data3(:,1), data3(:,6), 'LineWidth', 1.2, ...
                                                'color', [1 0.7 0])
% hold on
for i=2:O2.num_vibr_levels(1)-1
%     pause
    loglog(data3(:,1), data3(:,5+i), 'LineWidth', 1.2, ...
    'color', [1 0.7*(1-i/O2.num_vibr_levels(1)) i/O2.num_vibr_levels(1)])
% disp(i)
end
title('nO2(X, i)')
loglog(data3(:,1), data3(:,5+28), 'LineWidth', 1.2, ...
    'color', [0.1 0.8 0.1])
end
%%
function dddraw(x1, y1, x2, y2)
semilogx(x1, y1, 'LineWidth', 1.5)
hold on 
semilogx(x2, y2, 'LineWidth', 1.2)
legend('Tiago', 'Maksim', 'location', 'best')
xlim([0 x2(end)])
end