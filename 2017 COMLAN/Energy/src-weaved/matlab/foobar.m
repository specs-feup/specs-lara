function foobar()
   matisse_energy_start_0 = getCurrentEnergy();
   foo();
   matisse_energy_end_0 = getCurrentEnergy();
   matisse_energy_0 = (matisse_energy_end_0 - matisse_energy_start_0) / 1000000;
   fprintf('Energy:%f\n', matisse_energy_0);
   matisse_energy_start_1 = getCurrentEnergy();
   other();
   matisse_energy_end_1 = getCurrentEnergy();
   matisse_energy_1 = (matisse_energy_end_1 - matisse_energy_start_1) / 1000000;
   fprintf('Energy:%f\n', matisse_energy_1);
end

function [total] = getCurrentEnergy()
   % Check if folder exists
   raplFolder = '/sys/class/powercap/intel-rapl/';
   if (exist(raplFolder) ~= 7)
      disp(strcat('RAPL folder is not available: ', raplFolder));
      return;
   end
   % Get RAPL folders for each package
   %listing = dir(raplFolder);
   %listing
   raplFolderListing = eval('ls(raplFolder)');
   folders = strsplit(raplFolderListing);
   filter = strncmp(folders, 'intel-rapl:', numel('intel-rapl:'));
   packageFolders = folders(filter);
   energyReadings = zeros(size(packageFolders));
   for i = 1:numel(packageFolders)
      energyFile = strcat(packageFolders{i}, '/energy_uj');
      energyFilepath = strcat(raplFolder, energyFile);
      [status, output] = system(cstrcat('cat ', energyFilepath));
      energyReadings(i) = str2num(output);
   end
   total = sum(energyReadings);
end

function [a] = foo()
   a = 0;
   for i = 0:5
      matisse_energy_start_2 = getCurrentEnergy();
      a = a + bar();
      matisse_energy_end_2 = getCurrentEnergy();
      matisse_energy_2 = (matisse_energy_end_2 - matisse_energy_start_2) / 1000000;
      fprintf('Energy:%f\n', matisse_energy_2);
   end
end

function [result] = bar()
   result = 1.0;
end

