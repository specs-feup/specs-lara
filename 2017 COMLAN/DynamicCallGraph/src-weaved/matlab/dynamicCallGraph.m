function dynamicCallGraph()
   global dcg_CallGraph
   fprintf('digraph CallGraph {\n');
   if (numel(dcg_CallGraph) >= 1 && dcg_CallGraph(1) ~= 0)
      fprintf('\tfoobar -> foo [label=\"%d\"];\n', dcg_CallGraph(1));
   end
   if (numel(dcg_CallGraph) >= 2 && dcg_CallGraph(2) ~= 0)
      fprintf('\tfoobar -> other [label=\"%d\"];\n', dcg_CallGraph(2));
   end
   if (numel(dcg_CallGraph) >= 3 && dcg_CallGraph(3) ~= 0)
      fprintf('\tfoo -> bar [label=\"%d\"];\n', dcg_CallGraph(3));
   end
   fprintf('}\n');
   % Clear graph
   dcg_CallGraph = [];
end

