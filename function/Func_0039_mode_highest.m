function mode_value = Func_0039_mode_highest(slice_data, edges)
    % edges are normally -1:0.1:1

    if isfinite(mean(slice_data))
        counts = histcounts(slice_data(:), edges);
        
        [~, idx] = max(counts);
        mode_interval = [edges(idx), edges(idx+1)];
        
        mode_values = slice_data(slice_data >= mode_interval(1) & slice_data < mode_interval(2));
        mode_value = max(mode_values);
        mode_value=round(mode_value,1)-0.05;
        
        if mode(slice_data)==1
            mode_value=1;
        elseif mode(slice_data)==-1
            mode_value=-1;
        end
    else
        mode_value=NaN ;
    end

end