function logger_test()
   log_file_0 = fopen('log.txt', 'at+');
   fprintf(log_file_0, 'Logging to a file\n');
   fprintf('Print double %f after foo\n', 2.0);
   foo();
   fprintf('Printing again\n');
   fprintf(log_file_0, 'Logging again to a file\n');
   fprintf('Print double %f after bar\n', 2.0);
   fclose(log_file_0);
end

function [a] = foo()
   log_file_1 = fopen('log.txt', 'at+');
   a = 0;
   for i = 0:5
      fprintf(log_file_1, 'Logging to a file\n');
      fprintf('Print double %f after foo\n', 2.0);
      a = a + bar();
      fprintf('Printing again\n');
      fprintf(log_file_1, 'Logging again to a file\n');
      fprintf('Print double %f after bar\n', 2.0);
   end
   fclose(log_file_1);
end

function [result] = bar()
   result = 1.0;
end

