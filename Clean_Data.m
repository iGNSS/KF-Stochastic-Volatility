% Data should be cleaned from NaNs and Inf
%  Input: Data_to_clean - prices (should be a vector)
%         datas         - the date ticks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data_to_clean,datas] = Clean_Data(Data_to_clean,datas)
 index_2reject = find(isfinite(Data_to_clean)==0);   % find Inf and NaN index
 
 Data_to_clean(index_2reject) = [];    % delete bad entries
 datas(index_2reject,:) = [];          % the same for data
 disp([' -- The number of rejected data entries (NaN/INF) is ', num2str( length(index_2reject),'%8.1f' )]);
end
