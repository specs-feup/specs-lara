#define NUM_ATOMS 1000


double Distance2(double a[3], double b[3])
{
  return (a[0] - b[0]) * (a[0] - b[0]) +
         (a[1] - b[1]) * (a[1] - b[1]) +
         (a[2] - b[2]) * (a[2] - b[2]);
}

double SumOfInternalDistances(double atoms[][3], int numAtoms)
{
  int i, j;

  double total_dist = 0.0;

  for (int i = 0; i < numAtoms; i++) {
    for (int j = i; j < numAtoms; j++) {
      total_dist += Distance2(atoms[i], atoms[j]);
    }
  }
  return total_dist;
}


int main() {
	
	double atoms[NUM_ATOMS][3];		
	
	SumOfInternalDistances(atoms, NUM_ATOMS);
}