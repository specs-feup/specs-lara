function [b] = libcal(a)
   % monitoring b
   b = a + 1;
   monitor_range_update(3, b);
end

