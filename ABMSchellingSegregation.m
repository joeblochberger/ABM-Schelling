% ABMSchellingSegregation.m
% Tested on the following builds:  MATLAB R2020b
%
% (C) 2021 Joseph Blochberger (jblochb2@jhu.edu).  All rights reserved.
% Feel free to share and use this script with associated files based on the
% following creative commons license: Attribution-NonCommercial-ShareAlike
% 4.0 International (CC BY-NC-SA 4.0).  For more information, see
% creativecommons.org/licenses/by-na-sa/4.0/
%
% Kindly cite as
%       Blochberger, J. 2021. Agent-Based Modeling of Schelling's
%       Segregation Model. https://github.com/joeblochberger/ABM-Schelling/ABMSchellingSegregation.m, GitHub. Retrieved Month Day, Year.
%
% Usage:
% This script demonstrates a realization of Schelling's segregation model [1].
% This model has fixed boundary conditions and updates sequentially.
% The default demo is based on default inputs from NetLogo v6.2 [2].
% Note: user will need to define where output results go in addition to
% defining simulation start number to rerun this script once a simulation
% is complete.
%
% Input:
%      start - user types in a starting value for the beginning of the
%      simulation
%
%      finish - user types in a finishing value for the beginning of the
%      simulation
%
%      ScrptPath - user defined directory where this MATLAB script lives
%
%      ResPath - user defined directory where figures are populated along
%      with .mat files of results
%
% References:
%      [1] Schelling, T., 1971, “Dynamic Models of Segregation,” Journal of
%      Mathematical Sociology, 1, pp. 143-186.
%
%      [2] Wilensky, U. 1999. NetLogo. http://ccl.northwestern.edu/netlogo/
%      Center for Connected Learning and Computer-Based Modeling,
%      Northwestern University. Evanston, IL.
%
clc; close all; clear iter;

%% Setup
if exist('start','var') == 0
    clc; close all; clear all;
    ScrptPath = uigetdir('','Where is this script located?');
    ResPath = uigetdir('','Where would you like your results?');
    start = input('Enter the starting simulation number (e.g. 1): ');
    finish = input('Enter the final number of simulations you would like to run (e.g. 10): ');
    spyflag = input('Suppress visualizing unhappy agents? (Enter 1 to suppress, 0 to visualize): ');
    cd(ResPath);
    
    %% Inputs to ABM
    N = 25; % for a square lattice, length along one direction
    density = 95/100; % What percentage of the lattice do agents occupy?
    des_pct_sim = 50/100; % desired percent similar in neighborhood
    msz = 8; % plot markersize for unhappy agents
    
    % % % 1st order Moore neighborhood
    SurroundingNeighbors = [1,1,1;
        1,0,1;
        1,1,1];
    
elseif exist('start','var') == 1
    cd(ResPath);
    if SimulationNum >= finish
        clear all; close all;
        error('Simulations are finished.  Please check your results directory.');
    end
end

%% Simulation forloop
for SimulationNum = start:finish % Number of simulations to loop through
    % Initial population
    LatticeNetwork = zeros(N); % set up lattice
    agents = randi([0 1],size(LatticeNetwork,1),size(LatticeNetwork,2)); % populate agents on lattice network
    agents(agents==0) = 2; % differentiate some of them from the empty values
    empty_idx = randi([1 size(LatticeNetwork,1)*size(LatticeNetwork,2)], 1, round(size(LatticeNetwork,1).*size(LatticeNetwork,2)*(1-density))); %account for the predescribed density
    agents(empty_idx) = 0; % designate the empty values
    popagents = nonzeros(agents); total_agents=length(popagents);% total number of populated agents
    group_one = sum(popagents==1); % group one check
    group_two = sum(popagents==2); % group two check
    t = 0; % iteration zero
    
    % initial grid
    figure
    imagesc(agents);
    mymap = [1 1 1;1 0 0;0 1 0];colormap(mymap);
    labels={'Empty','Red','Green'};
    lcolorbar(labels,'fontweight','bold');
    hold on; xticks([]); yticks([]);
    
    %% Find neighbors in lattice using 2D convolution
    M = zeros(size(agents));
    pct_sim_mtx = zeros(size(agents));
    for x = 1:size(agents,1)
        for y = 1:size(agents,2)
            M(x,y)=1;
            nbrs{x,y}=agents(conv2(M,SurroundingNeighbors,'same')>0);
            if agents(x,y) == 1 % Red
                pct_sim_mtx(x,y) = sum(nbrs{x,y}==1)/length(nonzeros(nbrs{x,y}));
            elseif agents(x,y) == 2 % Green
                pct_sim_mtx(x,y) = sum(nbrs{x,y}==2)/length(nonzeros(nbrs{x,y}));
            elseif agents(x,y) == 0 % space is empty
                pct_sim_mtx(x,y) = 1;
            end
            M(:) = 0;
        end
    end
    
    % find the dissatisfied agents
    dis_mtx = pct_sim_mtx<des_pct_sim;
    num_unhappy = numel(nonzeros(dis_mtx));
    if spyflag == 0
        spy(sparse(dis_mtx),'kx',msz)
    end
    axes.YDir = 'reverse';
    pct_dis = sum(dis_mtx,'all')/length(popagents);
    fprintf(['Iteration = ',num2str(t),'\n'])
    fprintf(['Number Unhappy = ',num2str(num_unhappy),'\n'])
    fprintf(['Percentage Unhappy = ',num2str(pct_dis*100),'\n\n'])
    title('Iteration = 0')
    axis tight
    print([num2str(SimulationNum),'_SchellingSim_Iter',num2str(t)], '-dpng')
    close
    
    % document changes
    AgentsOverTime(1,:,:) = agents;
    num_unhappy(1) = num_unhappy;
    t=1; iter(1)=0;
    
    % find where unhappy agents can move to
    while pct_dis ~=0
        t = t+1;
        iter(t) = t-1;
        for x = 1:size(agents,1)
            for y = 1:size(agents,2)
                current_agent = agents(x,y);
                dis_agents = agents.*dis_mtx;
                dis_map = find(dis_agents~=0);
                empty_map = find(agents==0);
                
                if current_agent == dis_agents(x,y)
                    idx_2_mv = randi([1 length(empty_map)],1); % index to move arbitrarily to empty space
                    agents(empty_map(idx_2_mv)) = current_agent; % agent moves
                    agents(x,y) = 0; % to account for the agent moving
                    
                    % Find neighbors in lattice using 2D convolution after
                    % moving
                    M(x,y)=1;
                    nbrs{x,y}=agents(conv2(M,SurroundingNeighbors,'same')>0);
                    if agents(x,y)==1
                        pct_sim_mtx(x,y) = sum(nbrs{x,y}==1)/length(nonzeros(nbrs{x,y}));
                    elseif agents(x,y)==2
                        pct_sim_mtx(x,y) = sum(nbrs{x,y}==2)/length(nonzeros(nbrs{x,y}));
                    else
                        pct_sim_mtx(x,y) = 1;
                    end
                    M(:) = 0;
                else
                    continue
                end
            end
        end
        
        figure
        imagesc(agents);
        colormap(mymap);
        lcolorbar(labels,'fontweight','bold');
        hold on; xticks([]); yticks([]);
        M = zeros(size(agents));
        pct_sim_mtx = zeros(size(agents));
        for x = 1:size(agents,1)
            for y = 1:size(agents,2)
                M(x,y)=1;
                nbrs{x,y}=agents(conv2(M,SurroundingNeighbors,'same')>0);
                if agents(x,y)==1
                    pct_sim_mtx(x,y) = sum(nbrs{x,y}==1)/length(nonzeros(nbrs{x,y}));
                elseif agents(x,y)==2
                    pct_sim_mtx(x,y) = sum(nbrs{x,y}==2)/length(nonzeros(nbrs{x,y}));
                else
                    pct_sim_mtx(x,y) = 1;
                end
                M(:) = 0;
            end
        end
        
        % find the dissatisfied agents
        AgentsOverTime(t,:,:) = agents;
        dis_mtx = pct_sim_mtx<des_pct_sim;
        num_unhappy(t) = numel(nonzeros(dis_mtx));
        if spyflag == 0
            spy(sparse(dis_mtx),'kx',msz)
        end
        axes.YDir = 'reverse';
        pct_dis = sum(dis_mtx,'all')/length(popagents);
        fprintf(['Iteration = ',num2str(t-1),'\n'])
        fprintf(['Number Unhappy = ',num2str(num_unhappy(t)),'\n'])
        fprintf(['Percentage Unhappy = ',num2str(pct_dis*100),'\n\n'])
        title(['Iteration = ',num2str(t-1)])
        axis tight
        print([num2str(SimulationNum),'_SchellingSim_Iter',num2str(t-1)], '-dpng')
        close
        pause(0.5)
        
        if t==1000
            fprintf('No converged solution for unhappy agents \n\n')
            figure
            plot(iter,num_unhappy,'k','linewidth',2);
            title('Evolution of unhappiness')
            ylabel('Number of unhappy agents'); xlabel('Number of iterations')
            axis tight; grid on;
            print([num2str(SimulationNum),'_EvolutionOfUnhappiness'], '-dpng')
            close
            save([num2str(SimulationNum),'_SchellingSim.mat'],'AgentsOverTime','total_agents','iter','num_unhappy')
            fprintf('End of simulation \n\n')
            close all;
            clear agents axes current_agent dis_agents dis_map AgentsOverTime;
            clear dis_mtx empty_idx empty_map group_one group_two idx_2_mv;
            clear iter labels LatticeNetwork M mymap n_changes nbrs num_unhappy;
            clear numNeighbors pct_dis pct_sim_mtx popagents x y t;
            % Next simulation
            start = SimulationNum +1;
            if SimulationNum == finish
                clear all; close all;
                error('Simulations are finished.  Please check your results directory.');
                return
            elseif SimulationNum < finish
                run([ScrptPath,'\ABMSchellingSegregation.m'])
            end
        end
    end
    
    %%
    figure
    plot(iter,num_unhappy,'k','linewidth',2);
    title('Evolution of unhappiness')
    ylabel('Number of unhappy agents'); xlabel('Number of iterations')
    axis tight; grid on;
    print([num2str(SimulationNum),'_EvolutionOfUnhappiness'], '-dpng')
    close
    save([num2str(SimulationNum),'_SchellingSim.mat'],'AgentsOverTime','total_agents','iter','num_unhappy')
    fprintf('End of simulation \n\n')
    close all;
    clear agents axes current_agent dis_agents dis_map AgentsOverTime;
    clear dis_mtx empty_idx empty_map group_one group_two idx_2_mv;
    clear iter labels LatticeNetwork M mymap n_changes nbrs num_unhappy;
    clear numNeighbors pct_dis pct_sim_mtx popagents t x y;
    clc
end

% Next simulation
start = SimulationNum +1;
if SimulationNum == finish
    clear all; close all;
    error('Simulations are finished.  Please check your results directory.');
    return
elseif SimulationNum < finish
    run([ScrptPath,'\ABMSchellingSegregation.m'])
end
