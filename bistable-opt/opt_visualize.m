function stop = opt_visualize(x, optimValues, state)

    close all;

    stop = false;
    
    keypoints = opt_calculate(x);
    [rA1, rB1, rC1, rD1, rE1, rF1, rA2, rB2, rC2, rD2, rE2, rF2] = unpack_keypoints(keypoints);

    % Create figure & plotting settings
    fig = figure ();
    drawnow;
    MP = get(0, 'MonitorPositions');
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1:2) + MP(2, 1:2), pos(3:4)]);
    hold on;
    margin = 5;
    xArray = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
    yArray = horzcat(keypoints(2,:), keypoints(4,:)); % y-coordinates of all points of all joints
    xEnv = max(xArray) - min(xArray);
    yEnv = max(yArray) - min(yArray);

    % check which envelope dimension is larger
    if xEnv > yEnv
        greater = xEnv;
    else
        greater = yEnv;
    end

    % set upper limit using greater envelope dim to maintain equal axes
    xlim([min(xArray) - margin, min(xArray) + greater + margin])
    ylim([min(yArray) - margin, min(yArray) + greater + margin])
    
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_DF = plot([0],[0],'LineWidth',2);
    
    % Plot vectors between keypoints as lines
    set(h_OB, 'XData', [0, rB1(1)]);
    set(h_OB, 'YData', [0, rB1(2)]);
    
    set(h_AC, 'XData', [rA1(1), rC1(1)]);
    set(h_AC, 'YData', [rA1(2), rC1(2)]);
    
    set(h_BD, 'XData', [rB1(1), rD1(1)]);
    set(h_BD, 'YData', [rB1(2), rD1(2)]);
    
    set(h_CD, 'XData', [rC1(1), rD1(1)]);
    set(h_CD, 'YData', [rC1(2), rD1(2)]);
    
    set(h_CE, 'XData', [rC1(1), rE1(1)]);
    set(h_CE, 'YData', [rC1(2), rE1(2)]);
    
    set(h_DF, 'XData', [rD1(1), rF1(1)]);
    set(h_DF, 'YData', [rD1(2), rF1(2)]);

    % create 2nd set of axes so 2nd solution will plot with same color order
    a2 = axes;
    set(a2,'Visible','off') % hide axes background
    axes(a2)
%     xlim([-axlim, axlim])
%     ylim([-axlim, axlim])
    xlim([min(xArray) - margin, min(xArray) + greater + margin])
    ylim([min(yArray) - margin, min(yArray) + greater + margin])
    hold on;
    
    % Create plotting handles
    h_OB2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_AC2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_BD2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_CD2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_CE2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_DF2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    
    % Plot vectors between keypoints as lines
    set(h_OB2, 'XData', [0, rB2(1)]);
    set(h_OB2, 'YData', [0, rB2(2)]);
    
    set(h_AC2, 'XData', [rA2(1), rC2(1)]);
    set(h_AC2, 'YData', [rA2(2), rC2(2)]);
    
    set(h_BD2, 'XData', [rB2(1), rD2(1)]);
    set(h_BD2, 'YData', [rB2(2), rD2(2)]);
    
    set(h_CD2, 'XData', [rC2(1), rD2(1)]);
    set(h_CD2, 'YData', [rC2(2), rD2(2)]);
    
    set(h_CE2, 'XData', [rC2(1), rE2(1)]);
    set(h_CE2, 'YData', [rC2(2), rE2(2)]);
    
    set(h_DF2, 'XData', [rD2(1), rF2(1)]);
    set(h_DF2, 'YData', [rD2(2), rF2(2)]);

    % Label keypoints
    xtxt = [0, rA1(1), rB1(1), rC1(1), rD1(1), rE1(1), rF1(1), ...
        rC2(1), rE2(1), rF2(1)];
    ytxt = [0, rA1(2), rB1(2), rC1(2), rD1(2), rE1(2), rF1(2), ...
        rC2(2), rE2(2), rF2(2)];
    str = {'O', 'A', 'B', 'C1', 'D', 'E1', 'F1', 'C2', 'E2', 'F2'};
    text(xtxt,ytxt,str);

    hold off;

    % Save figure
    dir = strcat('images/', string(datetime('today', 'Format', 'yyyy-MM-dd')), '/');
    if ~exist(dir, 'dir')
        mkdir(dir);
    end
    iter = optimValues.iteration;
    saveas(fig, strcat(dir, num2str(iter), '.png'))
    close(fig)

    %{
    %% individual figures
    % plot first figure alone
    figure(2);
    hold on;
    axlim = 50;
%     xlim([-axlim, axlim])
%     ylim([-axlim, axlim])
    xlim([min(xArray) - margin, min(xArray) + greater + margin])
    ylim([min(yArray) - margin, min(yArray) + greater + margin])
    
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_DF = plot([0],[0],'LineWidth',2);
    
    % Plot vectors between keypoints as lines
    set(h_OB, 'XData', [0, rB1(1)]);
    set(h_OB, 'YData', [0, rB1(2)]);
    
    set(h_AC, 'XData', [rA1(1), rC1(1)]);
    set(h_AC, 'YData', [rA1(2), rC1(2)]);
    
    set(h_BD, 'XData', [rB1(1), rD1(1)]);
    set(h_BD, 'YData', [rB1(2), rD1(2)]);
    
    set(h_CD, 'XData', [rC1(1), rD1(1)]);
    set(h_CD, 'YData', [rC1(2), rD1(2)]);
    
    set(h_CE, 'XData', [rC1(1), rE1(1)]);
    set(h_CE, 'YData', [rC1(2), rE1(2)]);
    
    set(h_DF, 'XData', [rD1(1), rF1(1)]);
    set(h_DF, 'YData', [rD1(2), rF1(2)]);

    xtxt = [0, rA1(1), rB1(1), rC1(1), rD1(1), rE1(1), rF1(1)];
    ytxt = [0, rA1(2), rB1(2), rC1(2), rD1(2), rE1(2), rF1(2)];
    str = {'O', 'A', 'B', 'C1', 'D1', 'E1', 'F1'};
    text(xtxt,ytxt,str);

    hold off;


    % plot second solution alone
    figure(3)
%     xlim([-axlim, axlim])
%     ylim([-axlim, axlim])
    xlim([min(xArray) - margin, min(xArray) + greater + margin])
    ylim([min(yArray) - margin, min(yArray) + greater + margin])
    hold on;
    
    % Create plotting handles
    h_OB2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_AC2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_BD2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_CD2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_CE2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    h_DF2 = plot([0],[0],'LineWidth',2,'LineStyle','--');
    
    % Plot vectors between keypoints as lines
    set(h_OB2, 'XData', [0, rB2(1)]);
    set(h_OB2, 'YData', [0, rB2(2)]);
    
    set(h_AC2, 'XData', [rA2(1), rC2(1)]);
    set(h_AC2, 'YData', [rA2(2), rC2(2)]);
    
    set(h_BD2, 'XData', [rB2(1), rD2(1)]);
    set(h_BD2, 'YData', [rB2(2), rD2(2)]);
    
    set(h_CD2, 'XData', [rC2(1), rD2(1)]);
    set(h_CD2, 'YData', [rC2(2), rD2(2)]);
    
    set(h_CE2, 'XData', [rC2(1), rE2(1)]);
    set(h_CE2, 'YData', [rC2(2), rE2(2)]);
    
    set(h_DF2, 'XData', [rD2(1), rF2(1)]);
    set(h_DF2, 'YData', [rD2(2), rF2(2)]);

    xtxt = [0, rA2(1), rB2(1), rC2(1), rD2(1), rE2(1), rF2(1)];
    ytxt = [0, rA2(2), rB2(2), rC2(2), rD2(2), rE2(2), rF2(2)];
    str = {'O', 'A', 'B', 'C2', 'D2', 'E2', 'F2'};
    text(xtxt,ytxt,str);

    hold off;
    %}

end