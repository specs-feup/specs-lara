function [b] = libcal(a)
   % monitoring b
   b = a + 1;
   monitor1_range_update(3, b);
end

