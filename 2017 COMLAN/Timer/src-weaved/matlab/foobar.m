function foobar()
   tic;
   foo();
   matisse_time_0 = toc;
   fprintf('Time:%fus\n', matisse_time_0);
end

function [a] = foo()
   a = 0;
   for i = 0:5
      tic;
      a = a + bar();
      matisse_time_1 = toc;
      fprintf('Time:%fus\n', matisse_time_1);
   end
end

function [result] = bar()
   result = 1.0;
end

