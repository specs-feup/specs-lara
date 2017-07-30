function monitor1_range_update(index, value)
   global monitor1_range_min monitor1_range_max
   if (value < monitor1_range_min(index))
      monitor1_range_min(index) = value;
   end
   if (value > monitor1_range_max(index))
      monitor1_range_max(index) = value;
   end
end

