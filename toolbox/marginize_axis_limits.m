%
% Marginize axis limit
%
% Constrain axis limits to a given margin
%
% Inputs:
%   limits = [v_min, v_max]
%   margin = constraint margin about the limits midpoint
%
% Outputs:
%   marginized_limits = [v_min, v_max]
%

function marginized_limits = marginize_axis_limits(limits, margin)
    
    v_min = limits(1);
    v_max = limits(2);
    
    v_margin = 0.5 * (v_max - v_min);
    if v_margin < margin
        v_mid = v_min + v_margin;
        v_min = v_mid - margin;
        v_max = v_mid + margin;
    end
    
    marginized_limits = [v_min, v_max];
    
end
