function foobar()
   monitor_range_init();
   foo();
   monitor_range_print();
end

function [a] = foo()
   % monitoring a
   a = 0;
   monitor_range_update(1, a);
   for i = 0:5
      % monitoring a
      a = a + bar();
      monitor_range_update(1, a);
   end
   % monitoring a
   a = libcall(a);
   monitor_range_update(1, a);
end

function [result] = bar()
   % monitoring result
   result = 1.0;
   monitor_range_update(2, result);
end

