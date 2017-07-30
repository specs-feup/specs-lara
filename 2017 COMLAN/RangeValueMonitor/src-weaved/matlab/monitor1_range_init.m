function monitor1_range_init()
   global monitor1_range_min monitor1_range_max
   for i = 1:3
      monitor1_range_min(i) = Inf;
      monitor1_range_max(i) = -Inf;
   end
end

