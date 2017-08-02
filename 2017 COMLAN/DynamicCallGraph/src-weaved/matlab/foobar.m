function foobar()
   global dcg;
   foo();
   if (numel(dcg) < 1)
      dcg(1) = 1;
   else
      dcg(1) = dcg(1) + 1;
   end
   other();
   if (numel(dcg) < 2)
      dcg(2) = 1;
   else
      dcg(2) = dcg(2) + 1;
   end
   dynamicCallGraph
end

function [a] = foo()
   global dcg;
   a = 0;
   for i = 0:5
      a = a + bar();
      if (numel(dcg) < 3)
         dcg(3) = 1;
      else
         dcg(3) = dcg(3) + 1;
      end
   end
end

function [result] = bar()
   result = 1.0;
end

