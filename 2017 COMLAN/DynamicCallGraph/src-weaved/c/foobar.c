#include <stdio.h>
int dcg_CallGraph[2] = {0};

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      a += bar();
      dcg_CallGraph[0]++;
   }
   
   return a;
}


int main() {
   foo();
   dcg_CallGraph[1]++;
   printf("digraph CallGraph {\n");
   	if (dcg_CallGraph[ 0 ] != 0) {
   		printf("\tfoo -> bar [label=\"%d\"];\n", dcg_CallGraph[0]);
   	}
   	if (dcg_CallGraph[ 1 ] != 0) {
   		printf("\tmain -> foo [label=\"%d\"];\n", dcg_CallGraph[1]);
   	}
   printf("}\n");
}
