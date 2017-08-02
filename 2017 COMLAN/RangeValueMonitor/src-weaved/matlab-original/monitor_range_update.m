function monitor_range_update(index, value)
   global monitor_range_min monitor_range_max
   if (value < monitor_range_min(index))
      monitor_range_min(index) = value;
   end
   if (value > monitor_range_max(index))
      monitor_range_max(index) = value;
   end
end

