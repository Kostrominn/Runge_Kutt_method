#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const double alpha = 0.001;
const double EPS = 1e-12;

typedef struct
{
    double x0;
    double y0;
}dos;
double max(double x, double y);
void coeffs(double * c, double ** a, double * b);
double * make (double * data, int n);
double f1(double t, double x, double y);
double f2(double t, double x, double y);
dos runge_Kutt(double h, double y1, double y2, double t, double * c, double ** a, double * b);
double norm(dos tmp1, dos tmp2);
void solvef(int n, double l, double r, double y1, double y2, double * c, double ** a, double * b);
void solvec(double l, double r, double y1, double y2, double eps, double * c, double ** a, double * b);
double Method_Hord( double t1, double t2, dos tmp1, dos tmp2, double * c, double ** a, double * b);
double Method_Hord_x( double t1, double t2, dos tmp1, dos tmp2, double * c, double ** a, double * b);
//double period (double* y);

int main()
{
    dos tmp1, tmp2, tmp;
    double * b;
    double * c;
    double ** a;

    int n,i;
    double eps;
    double h;
    double l, r, x0, y0;
    double t, y, yn;
    double s;
    printf("Input n\n");
    scanf("%d", &n);
    printf("Input l and r\n");
    scanf("%lf %lf", &l, &r);
    printf("Input y_0\n");
    scanf("%lf %lf", &x0, &y0);
    h = (r - l)/n;
    b = (double *)malloc(sizeof(double)*7);
    c = (double *)malloc(sizeof(double)*7);
    a = (double **)malloc(sizeof(double *)*6);

    if(!a || !b || !c)
    {
        printf("Cannot allocate memory \n");
        return 1;
    }
    for(i = 0; i < 6; i++)
        a[i] = (double *)malloc(sizeof(double)*(i + 1));
    coeffs(c, a, b);
    //solveFixed(n, l, r, y1_0, y2_0, c, a, b);
    solvec (l, r, x0, y0, EPS, c, a, b);
    free(c);
    free(b);
    for(i = 0; i < 6; i++)
        free(a[i]);
    free(a);
    return 0;
}

double max(double x, double y)
{
    if ( x > y ) return x;
    	else return y;
    //return (x > y) ? x : y;
}

double min(double x, double y)
{
    if ( x > y ) return y;
    	else return x;
    //return (x > y) ? x : y;
}
void coeffs(double * c, double ** a, double * b)
{
    c[0] = 0;
    c[1] = 0.5;
    c[2] = 2.0/3.0;
    c[3] = 1.0/3.0;
    c[4] = 5.0/6.0;
    c[5] = 1.0/6.0;
    c[6] = 1.0;
    b[0] = 13.0/200.0;
    b[1] = 0.0;
    b[2] = 11.0/40.0;
    b[3] = 11.0/40.0;
    b[4] = 4.0/25.0;
    b[5] = 4.0/25.0;
    b[6] = 13.0/200.0;
    a[0][0] = 1.0/2.0;
    a[1][0] = 2.0/9.0;
    a[1][1] = 4.0/9.0;
    a[2][0] = 7.0/36.0;
    a[2][1] = 2.0/9.0;
    a[2][2] = -1.0/12.0;
    a[3][0] = -35.0/144.0;
    a[3][1] = -55.0/36.0;
    a[3][2] = 35.0/48.0;
    a[3][3] = 15.0/8.0;
    a[4][0] = -1.0/360.0;
    a[4][1] = -11.0/36.0;
    a[4][2] = -1.0/8.0;
    a[4][3] = 1.0/2.0;
    a[4][4] = 1.0/10.0;
    a[5][0] = -41.0/260.0;
    a[5][1] = 22.0/13.0;
    a[5][2] = 43.0/156.0;
    a[5][3] = -118.0/39.0;
    a[5][4] = 32.0/195.0;
    a[5][5] = 80.0/39.0;
}
double * make (double * data, int n)
{
    double * tmp = (double *)malloc(sizeof(double)*2*n);
    memcpy(tmp, data, n*sizeof(double));
    free(data);
    return tmp;
}
double f1(double t, double x, double y)
{
    return y;
}
double f2(double t, double x, double y) 
{
    return (cos(t) - x - alpha * pow(x,3));
    //return -x;
}
double checkSol(double t)
{
    return t*t;
}
dos runge_Kutt(double h, double y1, double y2, double t, double * c, double ** a, double * b)
{
    dos output;
    double fin_1, fin_2;

    double k11, k21, k31, k41, k51, k61, k71;
    double k12, k22, k32, k42, k52, k62, k72; 
    
    //printf("h = %lf, y1 = %lf, y2 = %lf, t = %lf\n", h, y1, y2, t);
    k11 = f1(t, y1, y2);
    //printf("k11 = %lf\n", k11);
    k12 = f2(t, y1, y2);
    //printf("k12 = %lf\n", k12);
    k21 = f1(t + c[1]*h, y1 + a[0][0]*h*k11, y2 + a[0][0]*h*k12);
    //printf("k21 = %lf\n", k21);
    k22 = f2(t + c[1]*h, y1 + a[0][0]*h*k11, y2 + a[0][0]*h*k12);
    //printf("k22 = %lf\n", k22);
    k31 = f1(t + c[2]*h, y1 + a[1][0]*h*k11 + a[1][1]*h*k21, y2 + a[1][0]*h*k12 + a[1][1]*h*k22);
    //printf("k31 = %lf\n", k31);
    k32 = f2(t + c[2]*h, y1 + a[1][0]*h*k11 + a[1][1]*h*k21, y2 + a[1][0]*h*k12 + a[1][1]*h*k22);
    //printf("k32 = %lf\n", k32);
    k41 = f1(t + c[3]*h, y1 + a[2][0]*h*k11 + a[2][1]*h*k21 + a[2][2]*h*k31, y2 + a[2][0]*h*k12 + a[2][1]*h*k22 + a[2][2]*h*k32);
    k42 = f2(t + c[3]*h, y1 + a[2][0]*h*k11 + a[2][1]*h*k21 + a[2][2]*h*k31, y2 + a[2][0]*h*k12 + a[2][1]*h*k22 + a[2][2]*h*k32);
    k51 = f1(t + c[4]*h, y1 + a[3][0]*h*k11 + a[3][1]*h*k21 + a[3][2]*h*k31 + a[3][3]*h*k41, y2 + a[3][0]*h*k12 + a[3][1]*h*k22 + a[3][2]*h*k32 + a[3][3]*h*k42);
    k52 = f2(t + c[4]*h, y1 + a[3][0]*h*k11 + a[3][1]*h*k21 + a[3][2]*h*k31 + a[3][3]*h*k41, y2 + a[3][0]*h*k12 + a[3][1]*h*k22 + a[3][2]*h*k32 + a[3][3]*h*k42);
    k61 = f1(t + c[5]*h, y1 + a[4][0]*h*k11 + a[4][1]*h*k21 + a[4][2]*h*k31 + a[4][3]*h*k41 + a[4][4]*h*k51,  y2 + a[4][0]*h*k12 + a[4][1]*h*k22 + a[4][2]*h*k32 + a[4][3]*h*k42 + a[4][4]*h*k52);
    k62 = f2(t + c[5]*h, y1 + a[4][0]*h*k11 + a[4][1]*h*k21 + a[4][2]*h*k31 + a[4][3]*h*k41 + a[4][4]*h*k51,  y2 + a[4][0]*h*k12 + a[4][1]*h*k22 + a[4][2]*h*k32 + a[4][3]*h*k42 + a[4][4]*h*k52);
    //printf("k62 = %lf\n", k62);
    //exit(1);
    k71 = f1(t + c[6]*h, y1 + a[5][0]*h*k11 + a[5][1]*h*k21 + a[5][2]*h*k31 + a[5][3]*h*k41 + a[5][4]*h*k51 + a[5][5]*h*k61, y2 + a[5][0]*h*k12+ a[5][1]*h*k22 + a[5][2]*h*k32 + a[5][3]*h*k42 + a[5][4]*h*k52 + a[5][5]*h*k62);
    k72 = f2(t + c[6]*h, y1 + a[5][0]*h*k11 + a[5][1]*h*k21 + a[5][2]*h*k31 + a[5][3]*h*k41 + a[5][4]*h*k51 + a[5][5]*h*k61, y2 + a[5][0]*h*k12+ a[5][1]*h*k22 + a[5][2]*h*k32 + a[5][3]*h*k42 + a[5][4]*h*k52 + a[5][5]*h*k62);
    fin_1 = y1 + h*(b[0]*k11 + b[1]*k21 + b[2]*k31 + b[3]*k41 + b[4]*k51 + b[5]*k61 + b[6]*k71);
    fin_2 = y2 + h*(b[0]*k12 + b[1]*k22 + b[2]*k32 + b[3]*k42 + b[4]*k52 + b[5]*k62 + b[6]*k72);
    output.x0 = fin_1;
    output.y0 = fin_2;
    return output;
}   
double norm(dos tmp1, dos tmp2)
{
    return max(fabs(tmp1.x0 - tmp2.x0), fabs(tmp1.y0 - tmp2.y0));
    //return sqrt(pow(tmp1.y1 - tmp2.y1,2) + pow(tmp2.y2 - tmp1.y2, 2));
}
void solvef(int n, double l, double r, double y1, double y2, double * c, double ** a, double * b)
{
    double h;
    double period_v;
    double time_precise_0 = 0;
    double time_precise_1 = 0;
    dos tmp;
    double l_start = l;
    int sch = 1;
    //double sign_1;
    //double sign_2;

    FILE * out = fopen("check.dat", "w");
    h = (r - l)/n;
    tmp.x0 = y1;
    tmp.y0 = y2;

    //sign_1 = tmp.x0/abs(tmp.x0);
    printf("solvef\n");
    
    fprintf(out, "%lf %lf \n", tmp.y0, tmp.x0);
    while(l < r)
    {
        tmp = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);
        printf("tmp.x0 = %.10lf\n", tmp.x0);
        //printf("x[%d] = %lf\n", sch, x[sch]);
        //printf("y[%d] = %lf\n", sch, y[sch]);
        //printf("length[%d] = %lf\n", sch, length[sch]);
        //НЕ ФАКТ, ЧТО ЗДЕСЬ ЧЁТКО
        sch++;
        //sign_2 = tmp.x0/abs(tmp.x0);
        //if ( sign_2 * sign_1 < 0)
            //t_2 = l;
        //sign_1 = tmp.x0/abs(tmp.x0);

        l += h;
        fprintf(out, "%.10lf %.10lf \n", tmp.y0, tmp.x0);
    }
    
    //printf("h = %lf\n, n = %d, l = %lf\n", h, n, l);
    fclose(out);
}
void solvec(double l, double r, double y1, double y2, double eps, double * c, double ** a, double * b)
{
    double l_start = l;
    dos tmp, tmp1, tmp2, tmp3, tmp1p, tmp2p;
    double t1;
    double t2;
    double period = -1;
    int i = 0;
    int j = 0;
    double time_precise_0 = 0;
    double time_precise_1 = 0;
    double x_precise_1 = -1;
    double x_precise_0 = 0;
    FILE * out = fopen("data_5.dat", "w");

    int sch = 1;
    double h_beg = 0.1;
    double h;
    
    tmp.x0 = y1;
    tmp.y0 = y2;
    //printf("IN SOLVEC\n");
    //printf("tmp.x0 = %lf, tmp.y0 = %lf\n", tmp.x0, tmp.y0);
    tmp1p = tmp;
    tmp2p = tmp;
    t1 = l;
    t2 = l;
    fprintf(out, "%lf %lf \n", tmp.y0, tmp.x0);
    while((r - l) > EPS)
    {
        
        h = h_beg;
        tmp1 = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);
	    tmp2 = runge_Kutt(h/2.0, tmp.x0, tmp.y0, l, c, a, b);
	    tmp2 = runge_Kutt(h/2.0, tmp2.x0, tmp2.y0, l + h/2.0, c, a, b);
	        
        while (norm(tmp1, tmp2)/63. > eps)
        {
            h /= 2.0;
            //h * fmin(1.5, fmax(0, pow((tol/err), 1./(6. + 1.))))
            //printf("h = %lf, l = %lf \n", h, l);
        	tmp1 = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);
	        tmp2 = runge_Kutt(h/2.0, tmp.x0, tmp.y0, l, c, a, b);
	        tmp2 = runge_Kutt(h/2.0, tmp2.x0, tmp2.y0, l + h/2.0, c, a, b);
	        
	        //printf("SOlVEC WHILE_1, tmp1.x0 = %lf, tmp2.x0 = %lf\n", tmp1.x0, tmp2.x0);
            //printf("norm(tmp1, tmp2)/63. = %.10lf\n", norm(tmp1, tmp2)/63.);
	        //printf("i = %d\n", i);
	        //printf("h = %lf\n\n", h);
	        //i++;
        }
        //i =0;
        while (norm(tmp1, tmp2)/63. < eps/4.)
        {
            h*=2;
            tmp1 = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);
            tmp2 = runge_Kutt(h/2.0, tmp.x0, tmp.y0, l, c, a, b);
            tmp2 = runge_Kutt(h/2.0, tmp2.x0, tmp2.y0, l + h/2.0, c, a, b);
            //printf("SOlVEC WHILE_2, tmp1.x0 = %lf, tmp2.x0 = %lf\n", tmp1.x0, tmp2.x0);
            //printf("In RARE case\n");
            //printf("norm(tmp1, tmp2)/63. = %.10lf\n", norm(tmp1, tmp2)/63.);
            //printf("i = %d\n", i);
            //printf("h = %lf\n\n", h);
            //i++;	
        }
        if (l + h > r)
        {
            h = r - l;
        }
//-----------------------------------------------------FINDING PERIOD-----------------------------------------------//
        t2 = t1;
        t1 = l+h;
        tmp2p = tmp1p;
        tmp1p = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);

        tmp1 = runge_Kutt(h, tmp.x0, tmp.y0, l, c, a, b);
        tmp = tmp1;
        
        //printf("tmp1p.y0 = %lf, tmp2p.y0 = %lf\n", tmp1p.y0, tmp2p.y0);
        //printf("tmp1p.y0 * tmp2p.y0 = %lf\n", tmp2p.y0 * tmp1p.y0);
        if (tmp2p.y0 * tmp1p.y0 < 0)
        {
            //printf("x[%d] = %lf\n", time, x[time] );
            //printf("1 CHANGE\n %lf\n", l+h);
            //printf("time_precise = %.10lf", Method_Hord(time0, c, a, b));
            j++;
            if ( j == 1)
            {
                time_precise_0 = Method_Hord(t2, t1, tmp2p, tmp1p, c, a, b);
                x_precise_0 = Method_Hord_x(t2, t1, tmp2p, tmp1p, c, a, b);
                printf("j = 1, period = %lf \n", time_precise_0);
            }
            if ( j == 2)
            {
                time_precise_1 = Method_Hord(t2, t1, tmp2p, tmp1p, c, a, b);
                x_precise_1 = Method_Hord_x(t2, t1, tmp2p, tmp1p, c, a, b);
                printf("j = 2, period = %lf \n", time_precise_1);
            }
            if ( j == 3)
            {
                time_precise_1 = Method_Hord(t2, t1, tmp2p, tmp1p, c, a, b);
                x_precise_1 = Method_Hord_x(t2, t1, tmp2p, tmp1p, c, a, b);
                printf("j = 3, period = %lf \n", time_precise_0);
                printf("x_precise_0 = %.10lf, x_precise_1 = %.10lf, x_precise_1 - x_precise_0 = %e\n",x_precise_0,
                x_precise_1, x_precise_1 - x_precise_0);

                if (fabs(x_precise_1 - x_precise_0) < 0.01)
                {
                    printf("I'm IN LAST IF\n");
                    period = time_precise_1 - time_precise_0;
                }
                else
                {
                    printf("I'm IN LAST ELSE\n");
                    j = 2;
                    //time_precise_0 = time_precise_1;
                    //x_precise_0 = x_precise_1;
                }
            }
        }
        
//------------------------------------------------------------------------------------------------------------------//
        //printf("x = %lf\n", tmp.x0);
        //printf("y = %lf\n", tmp.y0);
        //printf("l = %lf, r = %lf\n", l, r);
        l += h;
        h_beg = h;
        fprintf(out, "%.10lf %.10lf \n", tmp.y0, tmp.x0);
    }
    //fprintf(out, "%.10lf %.10lf \n", l, tmp.x0);
    //exit(1);
    //printf("h = %lf\n, n = %d, l = %lf\n", (r-l_start)/n, n, l_start);
    //printf("%lf\n", period_v * (r-l_start)/n + l_start);
    printf("x_precise_1 - x_precise_0 = %e\n", x_precise_1 - x_precise_0);

    printf("period = %.10lf\n", period);
    fclose(out);
}

double Method_Hord(double t1, double t2, dos tmp1, dos tmp2, double * c, double ** a, double * b)
{
    double y1, y2, x1, t;
    double i = 0;
    dos nextxy;
   
    
    printf(" length1 = %lf\n",  t1 );
    printf(" length2 = %lf\n",  t2 );
	y1 = tmp1.y0;
    y2 = tmp2.y0;

    x1 = tmp1.x0;
    //int i = 0;
    printf("%f %f %f %f \n", t1, t2, y1, y2);
    while (fabs(y1) > 1e-8 && fabs(y2) > 1e-8)
    {
       // i++;
    	printf(" IN WHILE\n");
        t = (-y1 * t2 + y2 * t1) / (y2 - y1);
        //printf("t = %f, real t = %f \n", t, -y1 * t2 + y2 * t1);
        //next(nextxy, t1, x1, y1, t - t1);
        nextxy = runge_Kutt(t-t1, x1, y1, t1, c, a, b);
        //printf("nextxy.x0 = %.10lf\n", nextxy.x0);
        if (nextxy.y0 * y1 > 0)
        {
            //printf("IN IF\n");
            //printf("t1 = %.10lf, t2 =  %.10lf, y1 = %.10lf, y2 = %.10lf \n\n", t1, t2, y1, y2);
            t1 = t;
            x1 = nextxy.x0;
            y1 = nextxy.y0;
        }
        else
        {
            //printf("IN ELSE\n");
            t2 = t;
            y2 = nextxy.y0;
            //printf("t1 = %.10lf, t2 = %.10lf, y1 = %.10lf, x2 = %.10lf\n\n", t1, t2, y1, y2);
        }
        
        //printf("%f %f %f \n", y1, y2, t);
        //printf(" IN END OF WHILE y1 = %.12lf, fabs(y2) = %.10lf, t= %.10lf\n\n", y1, fabs(y2), t);
    }
    //printf("y1 = %.10lf, fabs(y2) = %.10lf, t= %.10lf\n\n", y1, fabs(y2), t);
    //printf("!%d\n", i);
    return t;
}
double Method_Hord_x( double t1, double t2, dos tmp1, dos tmp2, double * c, double ** a, double * b)
{
    double y1, y2, x1, t;
    double i = 0;
    dos nextxy;
   
    
    printf(" length1 = %lf\n",  t1 );
    printf(" length2 = %lf\n",  t2 );
    y1 = tmp1.y0;
    y2 = tmp2.y0;

    x1 = tmp1.x0;
    //int i = 0;
    printf("%f %f %f %f \n", t1, t2, y1, y2);
    while (fabs(y1) > 1e-8 && fabs(y2) > 1e-8)
    {
       // i++;
        printf(" IN WHILE\n");
        t = (-y1 * t2 + y2 * t1) / (y2 - y1);
        //printf("t = %f, real t = %f \n", t, -y1 * t2 + y2 * t1);
        //next(nextxy, t1, x1, y1, t - t1);
        nextxy = runge_Kutt(t-t1, x1, y1, t1, c, a, b);
        //printf("nextxy.x0 = %.10lf\n", nextxy.x0);
        if (nextxy.y0 * y1 > 0)
        {
            //printf("IN IF\n");
            //printf("t1 = %.10lf, t2 =  %.10lf, y1 = %.10lf, y2 = %.10lf \n\n", t1, t2, y1, y2);
            t1 = t;
            x1 = nextxy.x0;
            y1 = nextxy.y0;
        }
        else
        {
            //printf("IN ELSE\n");
            t2 = t;
            y2 = nextxy.y0;
            //printf("t1 = %.10lf, t2 = %.10lf, y1 = %.10lf, x2 = %.10lf\n\n", t1, t2, y1, y2);
        }
        
        //printf("%f %f %f \n", y1, y2, t);
        //printf(" IN END OF WHILE y1 = %.12lf, fabs(y2) = %.10lf, t= %.10lf\n\n", y1, fabs(y2), t);
    }
    //printf("y1 = %.10lf, fabs(y2) = %.10lf, t= %.10lf\n\n", y1, fabs(y2), t);
    //printf("!%d\n", i);
    return x1;
}