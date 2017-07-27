function foobarfoo()
   foo(1);
end

function [a] = foo(seed)
   a = 0;
   for i = seed:5
      a = a + bar(i);
   end
end

function [result] = bar(a)
   if (a > 0)
      result = a;
      return;
   end
   result = foo(a + 1);
end

