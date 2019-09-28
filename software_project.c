#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "coeffs.h"

double** cross_pdt (double** n1, double** n2);
double dot_pdt (double** a, double** b);
double** line_dir_pt(double** m,double** A);
int is_coplanar (double** p1,double** m1,double** p2,double** m2);
double** line_intersect(double** p1,double** m1,double** p2,double** m2);
double** create_vec(double a, double b, double c);
double** scal_pdt(double** y, int m, int n, double k);
void printfunc(double** p,double** m, char option);
double **A1,**m1,**A2,**m2;

int main()
{
	A1 = create_vec(1,0,0);
	m1 = create_vec(-1,2,2);
	A2 = create_vec(0,0,0);
	m2 = create_vec(2,-1,2);
	double** A3_A = scal_pdt (create_vec(4,1,1), 3, 1, 2.0/9.0);
	double** m3_A = create_vec(2,2,-1);
	double** A3_B = scal_pdt (create_vec(2,-1,2), 3, 1, 2.0/9.0);
	double** m3_B = create_vec(2,2,-1);
	double** A3_C = scal_pdt (create_vec(2,0,1), 3, 1, 1.0/3.0);
	double** m3_C = create_vec(2,2,-1);
	double** A3_D = create_vec(0,0,0);
	double** m3_D = create_vec(2,2,-1);

	double** l1_p = line_dir_pt(m1,A1);
	double** l2_p = line_dir_pt(m2,A2);
	double** l3_A = line_dir_pt(m3_A,A3_A);
	double** l3_B = line_dir_pt(m3_B,A3_B);
	double** l3_C = line_dir_pt(m3_C,A3_C);
	double** l3_D = line_dir_pt(m3_D,A3_D);

	printfunc(A3_A, m3_A,'A');
	printfunc(A3_B, m3_B,'B');
	printfunc(A3_C, m3_C,'C');
	printfunc(A3_D, m3_D,'D');
	
	savetxt(l1_p,"l1_p.dat",3,10);
	savetxt(l2_p,"l2_p.dat",3,10);
	savetxt(l3_A,"l3_A.dat",3,10);
	savetxt(l3_B,"l3_B.dat",3,10);
	savetxt(l3_C,"l3_C.dat",3,10);
	savetxt(l3_D,"l3_D.dat",3,10);	
	
	return 0;
}

double** cross_pdt (double** n1, double** n2)
{
	double ** n3 = createMat(3,3);
	    n3[0][0] = 0           ;   n3[1][0] = -(n1[2][0]) ;     n3[2][0] = n1[1][0]    ;
        n3[0][1] = n1[2][0]    ;   n3[1][1] = 0           ;     n3[2][1] = -(n1[0][0]) ;
        n3[0][2] = -(n1[1][0]) ;   n3[1][2] = n1[0][0]    ;     n3[2][2] = 0           ;        
    return matmul(n3, n2, 3,3,1);  
}

double dot_pdt (double** a, double** b)
{
	double** c = transpose (a,3,1);
	return **matmul(c,b,1,3,1);
}

double** line_dir_pt(double** m,double** A)
{
	int len = 10;
	double** x_AB = createMat(3,len);
	double** lam = linspace(0,10,len);
	for(int i=0; i<len; i++)
	{
		x_AB[0][i] = *A[0] + *lam[i]*(*m[0]);
		x_AB[1][i] = *A[1] + *lam[i]*(*m[1]);
		x_AB[2][i] = *A[2] + *lam[i]*(*m[2]);
	}
	return x_AB;
}

int is_coplanar (double** p1,double** m1,double** p2,double** m2)
{
	return abs((dot_pdt(linalg_sub(p1,p2,3,1), cross_pdt(m1,m2)))) < 1e-16 ;
}

double** line_intersect(double** p1,double** m1,double** p2,double** m2)
{
	double s = dot_pdt(cross_pdt(linalg_sub(p1,p2,3,1),m2),cross_pdt(m1,m2)) / pow(linalg_norm(cross_pdt(m1,m2),3), 2);
	double** x = createMat(3,1);
	*x[0] = *p1[0] - *m1[0]*s;
	*x[1] = *p1[1] - *m1[1]*s;
	*x[2] = *p1[2] - *m1[2]*s;
	return x;
}

double** create_vec(double a, double b, double c)
{
	double** x = createMat(3,1);
	*x[0] = a;*x[1] = b;*x[2] = c;
	return x;
}

double** scal_pdt(double** y, int m, int n, double k)
{
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			y[i][j] = k*(y[i][j]);
		}
	}
	return y;
}

void printfunc(double** p,double** m, char option)
{
		int c1 = is_coplanar(A1, m1, p, m);
		int c2 = is_coplanar(A2, m2, p, m);
		int y = abs(linalg_norm(cross_pdt(m, cross_pdt(m1, m2)),3)) < 2e-16;
		
		if (!c1 || !c2 || !y)
			printf("option %c does not satisfy the condition\n", option); 
		else
			printf("option %c satisfies the given condition\n", option);
}





