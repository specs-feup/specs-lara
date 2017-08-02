function monitor1_print_ranges()
   log_file_0 = fopen('ranges_foobar.m.txt', 'at+');
   global monitor_range_min monitor_range_max
   fprintf(log_file_0, 'foo\n');
   fprintf(log_file_0, '\ta: {%f, %f}\n', monitor_range_min(1), monitor_range_max(1));
   fprintf(log_file_0, 'bar\n');
   fprintf(log_file_0, '\tresult: {%f, %f}\n', monitor_range_min(2), monitor_range_max(2));
   fprintf(log_file_0, 'libcal\n');
   fprintf(log_file_0, '\tb: {%f, %f}\n', monitor_range_min(3), monitor_range_max(3));
   % End of function
   fclose(log_file_0);
end

