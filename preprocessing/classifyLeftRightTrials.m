% filepath: classifyLeftRightTrials.m
function behavior_data = classifyLeftRightTrials(behavior_data, session)
% CLASSIFYLEFTRIGHTTRIALS - Classifies trials as left or right based on position data and saves results
%
% Usage:
%   behavior_data = classifyLeftRightTrials(behavior_data, session)
%
% Inputs:
%   behavior_data: Structure containing linearized theta maze tracking data
%   session: Session identifier for saving data
%
% Output:
%   behavior_data: Same structure with added left_right field in trials.alternation
%
% The function analyzes the last 20% of each trial to determine if the animal
% went to the left (1) or right (2) arm, then saves the updated structure.

    % Extract trial parameters
    trial_start = behavior_data.trials.alternation.start;
    trial_end = behavior_data.trials.alternation.end;
    nTrials = behavior_data.trials.alternation.nTrials;
    
    % Determine left vs right trials based on polar_theta values
    left_right = zeros(1, nTrials);
    for i = 1:nTrials
        % Find indices corresponding to trial timepoints
        idx_start = find(behavior_data.timestamps >= trial_start(i), 1, 'first');
        idx_end = find(behavior_data.timestamps <= trial_end(i), 1, 'last');
        
        if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
            fprintf('Warning: Trial %d has invalid timestamps - skipping\n', i);
            continue; % Skip problematic trials
        end
        
        % Sample the end of the trial (last 20% of timepoints)
        window_size = ceil(0.2 * (idx_end - idx_start + 1));
        window = (idx_end - window_size + 1) : idx_end;
        
        % Get polar_theta values at the end of the trial
        theta_values = behavior_data.position.polar_theta(window);
        mean_theta = mean(theta_values, 'omitnan');
        
        % Classify based on polar_theta
        % Left = 1, Right = 2
        if mean_theta < -5
            left_right(i) = 1; % Left trial
        elseif mean_theta > 5
            left_right(i) = 2; % Right trial
        else
            % If unclear, try using x coordinate as fallback
            mean_x = mean(behavior_data.position.x(window), 'omitnan');
            if mean_x < 0
                left_right(i) = 1; % Left trial
            else
                left_right(i) = 2; % Right trial
            end
        end
    end
    
    % Store results in trial structure
    behavior_data.trials.alternation.left_right = left_right;
    
    % Set state names
    behavior_data.stateNames.left_right = {'Left', 'Right'};
    
    % Report statistics
    left_count = sum(left_right == 1);
    right_count = sum(left_right == 2);
    fprintf('Trial classification: %d left trials, %d right trials, %d unclassified trials\n', ...
        left_count, right_count, nTrials - left_count - right_count);
    
    % Save the updated behavioral data
    fprintf('Saving behavioral data with trial classifications...\n');
    saveStruct(behavior_data, 'behavior', 'session', session);
end