function monitor_range_init()
   global monitor_range_min monitor_range_max
   for i = 1:3
      monitor_range_min(i) = Inf;
      monitor_range_max(i) = -Inf;
   end
end

