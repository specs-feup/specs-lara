int call_graph[2] = {0};

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      a += bar();
      call_graph[0]++;
   }
   
   return a;
}


int main() {
   foo();
   call_graph[1]++;
   printf("digraph call_graph {\n");
   
   	if (call_graph[ 0 ] != 0)
   		printf("\tfoo -> bar [label=\"%d\"];\n", call_graph[0]);
   				
   
   	if (call_graph[ 1 ] != 0)
   		printf("\tmain -> foo [label=\"%d\"];\n", call_graph[1]);
   				
   printf("}\n");
}
