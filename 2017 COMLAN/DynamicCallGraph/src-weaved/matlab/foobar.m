function foobar()
   global dcg_CallGraph;
   foo();
   if (numel(dcg_CallGraph) < 1)
      dcg_CallGraph(1) = 1;
   else
      dcg_CallGraph(1) = dcg_CallGraph(1) + 1;
   end
   other();
   if (numel(dcg_CallGraph) < 2)
      dcg_CallGraph(2) = 1;
   else
      dcg_CallGraph(2) = dcg_CallGraph(2) + 1;
   end
   dynamicCallGraph
end

function [a] = foo()
   global dcg_CallGraph;
   a = 0;
   for i = 0:5
      a = a + bar();
      if (numel(dcg_CallGraph) < 3)
         dcg_CallGraph(3) = 1;
      else
         dcg_CallGraph(3) = dcg_CallGraph(3) + 1;
      end
   end
end

function [result] = bar()
   result = 1.0;
end

