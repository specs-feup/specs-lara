function monitor1_print_ranges()
   global monitor_range_min monitor_range_max
   log_file = fopen('ranges_foobar.m.txt', 'at+');
   fprintf(log_file, 'foo\n');
   fprintf(log_file, '\tb: {%f, %f}\n', monitor_range_min(1), monitor_range_max(1));
   fprintf(log_file, 'bar\n');
   fprintf(log_file, '\tb: {%f, %f}\n', monitor_range_min(2), monitor_range_max(2));
   fprintf(log_file, 'libcal\n');
   fprintf(log_file, '\tb: {%f, %f}\n', monitor_range_min(3), monitor_range_max(3));
   fclose(log_file);
   % End of function
end

