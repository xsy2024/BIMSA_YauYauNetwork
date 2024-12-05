clear;clc;close all;
addpath('/Users/xsy/Documents/MATLAB/NonlinearFiltering-master/subprograms');
addpath('/Users/xsy/Documents/MATLAB/NonlinearFiltering-master/h');
%% Settings
Dim    = 3;
T      = 20;
dT     = 0.001;
dTau   = 5*dT;
dX     = 0.5;

% A1 = [-3.4554092 , -2.68897588,  7.60570567;
%     -1.97236709, -0.78826619,  3.68802256;
%     -0.82320704, -1.15490118,  1.8064858];
% A2 = [0.03093756,  0.07581195,-0.15806514;
%     0.01747811,  0.02475632,-0.06153576;
%     0.00687582,  0.03301364,-0.0320007];
% Ac = [0.35002403, -2.37092631, -0.19224301;
%     0.12984874, -1.23160753, -0.17537569;
%     0.07560191, -0.72405088,  0.00618002];
% A_0 = [33.41398988178064, 15.857458748688114, 11.827912779635522];
% f = @(x) (x*A1' + x.^2*A2'+ sin(x)*Ac'+ A_0);
% df = @(x) (diag(A1)'+2*x.*diag(A2)'+cos(x).*diag(Ac)');

f  = @(x)   [-0.00994309*sin(2*x(:,1))-0.00994309*cos(2*x(:,1))-0.00994309*sin(2*x(:,2))-0.01205387*sin(3*x(:,3))-0.01205387*cos(3*x(:,3)),...
            -0.04031238*sin(2*x(:,1))-0.04031238*cos(2*x(:,1))-0.04031238*sin(2*x(:,2))-0.03348778*sin(3*x(:,3))-0.03348778*cos(3*x(:,3)),...
            -0.03189308*sin(2*x(:,1))-0.03189308*cos(2*x(:,1))-0.03189308*sin(2*x(:,2))-0.03854743*sin(3*x(:,3))-0.03854743*cos(3*x(:,3))];
df = @(x)   [-0.00994309*2*cos(2*x(:,1))+0.00994309*2*sin(2*x(:,1)),...
            -0.04031238*2*cos(2*x(:,2)),...
             -0.03854743*3*cos(3*x(:,3))-0.03854743*3*sin(3*x(:,3))];
%  
% f = @(x) [-0.5*x(:,1)+0.1*x(:,2),...
%                0.1*x(:,3)-0.5*x(:,2), ...
%                -0.5*x(:,3)];
%  df = @(x) -0.5;
% f = @(x) 0.8*x.*(1-x./50);
% df = @(x) 0.016*(10-2*x);
% h = @(x) x;
 % f  = @(x)   [10*(x(:,2)-x(:,1)),...
 %             x(:,1).*(28-x(:,3))-x(:,2), ...
 %             x(:,1).*x(:,2)-8/3*x(:,3)];
 % df = @(x)   [-10,-1,-8/3];
 % f = @(x) [128.8245+0.266025*x(:,1).^2-0.33067886*x(:,2).^2-1.50997027*x(:,3).^2-0.00332551*x(:,1).^3+0.00870392*x(:,2).^3+...
 %         0.04036269*x(:,3).^3+0.42552532*cos(x(:,1))+1.90295568*cos(x(:,2))+1.85962133*cos(x(:,3)),...
 % 67.0436+0.14271644*x(:,1).^2-0.14707208*x(:,2).^2-0.84628994*x(:,3).^2-0.00179267*x(:,1).^3+0.00393847*x(:,2).^3+...
 %         0.02286477*x(:,3).^3+0.17546576*cos(x(:,1))+0.96474161*cos(x(:,2))+ 1.01792932*cos(x(:,3)),...
 % 44.7125+0.09310835*x(:,1).^2-0.11024018*x(:,2).^2-0.5288023*x(:,3).^2-0.00116232*x(:,1).^3+0.00288998*x(:,2).^3+...
 %         0.01408074*x(:,3).^3+0.14869103*cos(x(:,1))+0.52766047*cos(x(:,2))+0.75290533*cos(x(:,3))];
 % 
 % df = @(x) [2*0.266025*x(:,1)-3*0.00332551*x(:,1).^2-0.42552532*sin(x(:,1)),...
 %     -2*0.14707208*x(:,2)+3*0.00393847*x(:,2).^2-0.96474161*sin(x(:,2)),...
 %     -2*0.5288023*x(:,3)+3*0.01408074*x(:,3).^2-0.75290533*sin(x(:,3))];
%  f = @(x) [4.3465-0.581*x(:,1)-3.8668*x(:,2)+5.4308*x(:,3)+0.00507*x(:,1).^2+0.0707*x(:,2).^2-0.1142*x(:,3).^2,...
%      -16.583-1.244*x(:,1)-0.1176*x(:,2)+4.178*x(:,3)+0.010885*x(:,1).^2+0.00077785*x(:,2).^2-0.0795*x(:,3).^2,...
%      20.3121-0.5610*x(:,1)-2.3455*x(:,2)-0.3224*x(:,3)-0.00477105*x(:,1).^2+0.04497*x(:,2).^2+0.00272*x(:,3).^2
% 
% ];
%  df = @(x) [-0.581+2*0.00507*x(:,1),-0.1176+2*0.00077785*x(:,2),-0.3224+0.00272*2*x(:,3)];
 % f = @(x) [0.266025*x(:,1).^2-0.33067886*x(:,2).^2-1.50997027*x(:,3).^2+0.42552532*cos(x(:,1))+1.90295568*cos(x(:,2))+1.85962133*cos(x(:,3)),...
 % 0.14271644*x(:,1).^2-0.14707208*x(:,2).^2-0.84628994*x(:,3).^2+0.17546576*cos(x(:,1))+0.96474161*cos(x(:,2))+ 1.01792932*cos(x(:,3)),...
 % 0.09310835*x(:,1).^2-0.11024018*x(:,2).^2-0.5288023*x(:,3).^2+0.14869103*cos(x(:,1))+0.52766047*cos(x(:,2))+0.75290533*cos(x(:,3))];
 % 
 % df = @(x) [2*0.266025*x(:,1)-0.42552532*sin(x(:,1)),...
 %     -2*0.14707208*x(:,2)-0.96474161*sin(x(:,2)),...
 %     -2*0.5288023*x(:,3)-0.75290533*sin(x(:,3))];

 %g = @(x) [0.266025*x(:,1).^2-0.33067886*x(:,2).^2-1.50997027*x(:,3).^2-0.00332551*x(:,1).^3+0.00870392*x(:,2).^3+0.04036269*x(:,3).^3+0.42552532*cos(x(:,1))+1.90295568*cos(x(:,2))+1.85962133*cos(x(:,3))]

 % f  = @(x)   [-0.5*x(:,1)+cos(x(:,1))+0.1*x(:,2),...
 %              0.1*x(:,3)-0.5*x(:,2)+cos(x(:,2)), ...
 %              -0.5*x(:,3)+cos(x(:,3))];
 % df = @(x)   [-0.5-sin(x(:,1)),-0.5-sin(x(:,2)),-0.5-sin(x(:,3))];
 % 

 %  f  = @(x)   [-0.5*sin(3*x(:,1))+0.1*sin(2*x(:,2))+0.1*cos(x(:,1)),...
 %               -0.5*sin(3*x(:,2))+0.1*sin(2*x(:,3))+0.1*cos(x(:,2)),...
 %              -0.5*sin(3*x(:,3))+0.1*cos(x(:,3))];
 % df = @(x)   [-0.5*3*cos(3*x(:,1))-0.1*sin(x(:,1)),...
 %               -0.5*3*cos(3*x(:,2))-0.1*sin(x(:,2)),...
 %              -0.5*3*cos(x(:,3))-0.1*sin(x(:,3))];

 %h  = @(x)   x(:,1).^3+x(:,2).^3+x(:,3).^3;

 %h = @(x)  h_mat('/Users/xsy/Documents/MATLAB/NonlinearFiltering-master/h',x);
 %h  = @(x)   x(:,1)+x(:,2)+x(:,3);                                
 %h = @(x) x.^3 ;
 %h = @(x) [sqrt(1/3)*(x(:,1)+x(:,2)+x(:,3)),-2/3*x(:,1)+1/3*x(:,2)+1/3*x(:,3),(2/3*sqrt(1/3)-1/3)*x(:,1)+(-1/3*sqrt(1/3)-1/3)*x(:,2)+(-1/3*sqrt(1/3)+2/3)*x(:,3)];
 %h = @(x) [sqrt(1/3)*(x(:,1)+x(:,2)+x(:,3)),sqrt(1/2)*(x(:,1)-x(:,2)),sqrt(1/6)*(x(:,1)+x(:,2)-2*x(:,3))];
 %h = @(x) [sqrt(1/3)*(x(:,1)+x(:,2)+x(:,3)),sqrt(1/2)*(x(:,1)-x(:,2)),sqrt(1/6)*(x(:,1)+x(:,2)-2*x(:,3))];

%h = @(x) x*A(:,1:3)'+x.^2*A(:,4:6)'+cos(x)*A(:,7:9)'+B;

A = [1.11321727e+00, -1.95579459e+01,  2.34884903e+01;
    3.38276311e-01,  2.67941452e+00, -3.95948370e+00;
    6.26430811e-01, -1.04278484e+01,  1.20274497e+01];
B = orth(A);
Tab = B/A;
AA = [-6.24608628e-04,9.12986507e-01, -2.35293418e+00;
    -3.62400933e-03,-2.32331748e-01,  1.25506619e-01;
    -3.81168845e-03,5.16732420e-01, -1.19791538e+00];
BB = Tab*AA;
 %h = @(x) [4.588*x(:,1)-13.602*x(:,2)-10.783*x(:,3),-1.84*x(:,1)-5.89*x(:,2)-5.3660*x(:,3),1.840*x(:,1)-4.792*x(:,2)-4.474*x(:,3)];
 h = @(x) x*A' + x.^2*AA';
 %h = @(x) x.^3*B' + x.^6*BB' ;

  %h = @(x) x*A'+ x.^2*AA'+cos(x)*Ac'+A0 ;
  %h = @(x) x*B'+ x.^2*BB'+cos(x)*Bc'+ B0' ;

  
  %h = @(x)  (x*orth(A1)' + x.^2*A22'+sin(x)*Acc'+ A_00');
  %h = @(x) f(x)*Th';
  
 %h = @(x) (A*x' + AA*(x.^2)' +AAA*(x.^3)'+A0)';
 %h = @(x) x*B' +B0';
  %h = @(x) x*A' ;
 %h = @(x) (x+S0)*A + (x+S0).^2*AA +(x+S0).^3*AAA+ A0';
 %h = @(x) x*B' + x.^2*BB' +x.^3*BBB'+ B0';
 %h = @(x) [sqrt(1/3)*(x(:,1).^3+x(:,2).^3+x(:,3).^3),sqrt(1/2)*(x(:,1).^3-x(:,2).^3),sqrt(1/6)*(x(:,1).^3+x(:,2).^3-2*x(:,3).^3)];

 %% Initialize
nTau   = T/dTau;
Tau    = 0:dTau:T;
nT     = dTau/dT;

fprintf('========================================\n');
fprintf([' ' num2str(Dim) '-D Yau-Yau Method using QIEM with DST\n']);
fprintf('========================================\n');

%% Generate States
fprintf('Generating States ...');
tic
%v = readmatrix("Total_data.csv");
v_1 = readmatrix("0_Total_data_E.csv");
v_2 = readmatrix("0_Total_data_P.csv");
v_3 = readmatrix("0_Total_data_S.csv");
%y = (TAB*v)';
a = 0.:0.25:36;
a_new = 0.:dT:36;
v_new1 = interp1(a, v_1(64,:), a_new, 'linear');
v_new2 = interp1(a, v_2(64,:), a_new, 'linear');
v_new3 = interp1(a, v_3(64,:), a_new, 'linear');
y1 = vertcat(v_new1,v_new2,v_new3)';
[state, obser, s] = SimulateStateObser(T, dT, f, h, Dim, y1);
toc
y = obser(1:nT:end,:);
%v = readmatrix("Total_data.csv");
v_1 = readmatrix("0_Total_data_E.csv");
v_2 = readmatrix("0_Total_data_P.csv");
v_3 = readmatrix("0_Total_data_S.csv");
%y = (TAB*v)';
%y = y*TAB';
%y = readmatrix("Q_1 (2).csv");

rX = [min(min(state)), max(max(state))];
fprintf('Range(States) = [%f, %f].\n', rX(1), rX(2));

tstart = tic;
%% Construct Matrix
fprintf('Constructing Matrices ...\n');
x = (rX(1):dX:rX(2)+dX).';
[Lambda, B, x, n] = KolmogorovEW(Dim, dT, x, f, df, h);

%% Solve Kolmogorov Equations
sigma0 = exp( -10 * ( sum(x.^2, 2) ));%这里设置了初值分布
Iu     = zeros((nTau-1)*nT+1, Dim);
Idx    = 1;
U      = sigma0;
U      = U / sum(U);%对U进行归一化
Iu(Idx,:) = sum(U(:,ones(Dim,1)).*(h(x)), 1);%取期望
Eu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);
for jj = 1:nTau %这是一个大时间步（需要更新观测）
    Idx = Idx + 1;
	fprintf('Computing step %d ... ', Idx);tic
    if jj == 1
        tmp = y(jj,:);
    else
        tmp = y(jj,:) - y(jj-1,:);
    end
    %y = py.torch.tensor(x);
    %U = NormalizedExp( sum( h(py.torch.tensor(py.numpy.array(single(x)))).*tmp(ones(size(x,1),1), :) , 2) ) .* U;
	U = NormalizedExp( sum( 1*h(x).*tmp(ones(size(x,1),1), :) , 2) ) .* U;%这里输入观测值
    U = DST_Solver(Dim, Lambda, B, U, n);
    U = U / sum(U);
	Iu(Idx,:) = sum(U(:,ones(Dim,1)).*(h(x)), 1);
    Eu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);
    toc
	for ii = 2:nT%这是一个小时间步
        Idx = Idx + 1;
        fprintf('Computing step %d ... ', Idx);tic
        U = DST_Solver(Dim, Lambda, B, U, n); %这里才用DST方法离散化矩阵
        Iu(Idx,:) = sum(U(:,ones(Dim,1)).*(h(x)), 1);
        Eu(Idx,:) = sum(U(:,ones(Dim,1)).*x, 1);
        toc
	end
end
telapsed = toc(tstart);



%% Plot the Result
%PlotState(T, dT, state, Iu);

path = "./sim_exp_fig/";

Csvfilename = "Iu_nlss.csv";
%writematrix(Eu,path+Csvfilename);

Iu = readmatrix(path+Csvfilename);
Csvfilename = "state_nlss.csv";
state = readmatrix(path+Csvfilename);
TotalT = 0:0.001:20;

%Iu =  readmatrix("10_noise1.csv");
figure("Position",[100,100,1200,500]);
%subplot(1,2,1);
t = tiledlayout(2, 3, "Position", [0.1,0.2,0.8,0.7],'TileSpacing', 'none', 'Padding', 'none');  % 1 行 3 列，间隙为"无"
ax = gobjects(size(state, 2), 1);
yMin = inf;  % 设置初始为正无穷
yMax = -inf;  % 设置初始为负无穷
ylabel(t, 'State', 'FontSize', 42); 
xlabel(t, 'Time','FontSize',42);
%s = h(state);
%s1 = h(Eu);
%Error_RMS = sqrt(mean((sum((s - Iu).^2, 2))/Dim, 1));
%Error_RMS_1 = sqrt(mean((sum((s - s1).^2, 2))/Dim, 1));

%t.TileSpacing = [0,5];
% 第一次遍历，计算所有图的y轴范围
for ii = 1:size(state, 2)
    yMin = min(yMin, min(Iu(:, ii)));  % 更新 Iu 的最小值
    yMin = min(yMin, min(state(:, ii)));  % 更新 state 的最小值
    yMax = max(yMax, max(Iu(:, ii)));  % 更新 Iu 的最大值
    yMax = max(yMax, max(state(:, ii)));  % 更新 state 的最大值
end
for ii = 1:size(state,2)
    ax(ii) = nexttile;
    %ax(ii).TileSize = [ii,1];
    %figure(ii,"OuterPosition", [100,100,100,100]);
    %plot(TotalT, obser(:,ii), 'g-');hold on
    %plot(TotalT, s(:,ii), 'g-');hold on
    %plot(TotalT, s1(:,ii), 'b-');hold on
    %plot(TotalT, y1(1:4001,ii), 'k-'); hold on
    plot(TotalT, Iu(:,ii), 'r-',"LineWidth",2); hold on
    %plot(TotalT, Eu(:,ii), 'k-'); hold on
    plot(TotalT, state(:,ii), 'k-',"LineWidth",2);
    title("Variable "+ii,'FontSize', 36);
    %ylabel('state','FontSize',16);
    %legend('Estimates','States');
    ylim([yMin, yMax+0.5]);
    if ii > 1
        ax(ii).YColor = 'k';  % 隐藏 y 轴的颜色和刻度
        ax(ii).YTickLabel = [];
        ax(ii).YTick = ax(1).YTick; 
        ax(ii).YLabel = [];
    end
    if ii==1
        ax(ii).YTick = ax(1).YTick;
        disp(ax(1).YTickLabel(2:end));
        ax(ii).YTickLabel(1) = {""} ;
    end
    ax(ii).XTick = ax(1).XTick;
    ax(ii).XTickLabel = [];
    ax(ii).XLabel = [];
    set(ax(ii), "FontSize", 32, "LineWidth",2);
    if ii == 3
        annotation("textbox",[.85 .85 .1 .2],...
            "string", "Non-orth","EdgeColor","none","Rotation",270,'FontSize',32);
        %text(0.9, 0.5, "A", 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    end
       % if ii == 1
    %     ylabel('State', 'FontSize', 24);  % 显示 y 轴标签
    % end

    % if ii == 2
    %     xlabel('Time','FontSize',24);
    % end
end

%xlabel('time','FontSize',16);
linkaxes(ax, 'x');
Csvfilename = "IU_nlss_orth.csv";
%writematrix(Iu,path+Csvfilename);
% Csvfilename = "state_sindy.csv";
% writematrix(state,path+Csvfilename);

Iu = readmatrix(path+Csvfilename);
Csvfilename = "state_nlss_orth.csv";
state = readmatrix(path+Csvfilename);
yMin = inf;  % 设置初始为正无穷
yMax = -inf;  % 设置初始为负无穷

% 第一次遍历，计算所有图的y轴范围
for ii = 1:size(state, 2)
    yMin = min(yMin, min(Iu(:, ii)));  % 更新 Iu 的最小值
    yMin = min(yMin, min(state(:, ii)));  % 更新 state 的最小值
    yMax = max(yMax, max(Iu(:, ii)));  % 更新 Iu 的最大值
    yMax = max(yMax, max(state(:, ii)));  % 更新 state 的最大值
end
for ii = 1:size(state,2)
    ax(ii) = nexttile;
    %ax(ii).TileSize = [ii,1];
    %figure(ii,"OuterPosition", [100,100,100,100]);
    %plot(TotalT, obser(:,ii), 'g-');hold on
    %plot(TotalT, y1(1:4001,ii), 'k-'); hold on
    plot(TotalT, Iu(:,ii), 'r-',"LineWidth",2); hold on
    plot(TotalT, state(:,ii), 'k-',"LineWidth",2);
    %title("Variable "+ii,'FontSize', 24);
    %ylabel('state','FontSize',16);
    %legend('Estimates','States');
    ylim([yMin, yMax+0.5]);
    if ii > 1
        ax(ii).YColor = 'k';  % 隐藏 y 轴的颜色和刻度
        ax(ii).YTickLabel = [];
        ax(ii).YTick = ax(1).YTick; 
        ax(ii).YLabel = [];
    end
    
    if ii < 3
    ax(ii).XTick = ax(ii).XTick(1:end-1);
    else
    ax(ii).XTick = ax(3).XTick;
    end
    set(ax(ii), "FontSize", 32, "LineWidth",2);
    if ii == 3
        annotation("textbox",[.85 .45 .1 .2],...
            "string", "Orth","EdgeColor","none","Rotation",270,'FontSize',32);
    end
    % if ii == 1
    %     ylabel('State', 'FontSize', 24);  % 显示 y 轴标签
    % end
    % 
    % if ii == 2
    %     xlabel('Time','FontSize',24);
    % end
end

%xlabel('time','FontSize',16);
linkaxes(ax, 'x');

%Error_RMS = sqrt(mean((sum((s - Iu).^2, 2))/Dim, 1));
%Error_RMS = sqrt(mean((sum((s - s1).^2, 2))/Dim, 1));
%Error_M   = mean(sqrt(sum((state - Iu).^2, 2)/Dim), 1);
%Error_RMS = sqrt(mean((sum((y1(1:40001,:) - Iu).^2, 2))/Dim, 1));
%Error_M   = mean(sqrt(sum((y1(1:40001,:) - Iu).^2, 2)/Dim), 1);

%Csvfilename = "IU.csv";
%writematrix(Iu,Csvfilename);
 % path = "./sim_exp_fig/";
 % Csvfilename = "IU_lss_1.csv";
 % writematrix(Iu,path+Csvfilename);
 % Csvfilename = "state_lss_1.csv";
 % writematrix(state,path+Csvfilename);



fprintf('=============================================================\n');
fprintf(['The computation of ' num2str(Dim) '-D Yau-Yau Method is finished.\n']);
fprintf('-------------------------------------------------------------\n');
fprintf(['Terminal Time         : ' num2str(T) ' \n']);
fprintf(['Space Range           : [' num2str(rX(1)) ', ' num2str(rX(2)) '] \n']);
fprintf(['Size of Time Steps    : ' num2str(dT) ' \n']);
fprintf(['Size of Space Steps   : ' num2str(dX) ' \n']);
fprintf(['Root-Mean-Square Error: ' num2str(Error_RMS) ' \n']);
fprintf(['Mean Error            : ' num2str(Error_M) ' \n']);
fprintf(['Time Costs            : ' num2str(telapsed) ' seconds. \n']);
fprintf('=============================================================\n');
