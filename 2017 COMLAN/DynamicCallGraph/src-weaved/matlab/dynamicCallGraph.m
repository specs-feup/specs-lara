function dynamicCallGraph()
   global dcg
   fprintf('digraph dcg {\n');
   if (numel(dcg) >= 1 && dcg(1) ~= 0)
      fprintf('\tfoobar -> foo [label=\"%d\"];\n', dcg(1));
   end
   if (numel(dcg) >= 2 && dcg(2) ~= 0)
      fprintf('\tfoobar -> other [label=\"%d\"];\n', dcg(2));
   end
   if (numel(dcg) >= 3 && dcg(3) ~= 0)
      fprintf('\tfoo -> bar [label=\"%d\"];\n', dcg(3));
   end
   fprintf('}\n');
   % Clear graph
   dcg = [];
end

