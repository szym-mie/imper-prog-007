#include <stdio.h>
#include <math.h>

#define E 

#define TEST 1   // 1 - dla testowania, 0 - dla sprawdzarki

#define RECURS_LEVEL_MAX  10
#define N_MAX             10

/////////////////////////////////////////////////////////////////////
// 7.1 Calki jednokrotne

// Przyklady podcalkowych funkcji jednej zmiennej

typedef double (*Func1vFp)(double);  // Definiowanie nazwy dla typu wskaznika do funkcji jedej zmiennej 

double f_poly(double x) {  //    polynomial  a[0]+a[1]x+ ... + a[n]x^n
	return 2 * pow(x, 5) - 4 * pow(x, 4) + 3.5 * pow(x, 2) + 1.35 * x - 6.25;
}

double f_rat(double x) {  // Przyklad funkcji podcalkowej Proba z funkcja f(x) = 1/((x-0.5)*(x-0.5)+0.01) w przedziale od 0 do 3.
	return 1 / (pow(x-0.5, 2) + 0.01);
}

double f_exp(double x) {  // Przyklad funkcji podcalkowej
	return 2 * x * exp(-1.5 * x) - 1;
}

double f_trig(double x) {  // Przyklad funkcji podcalkowej
	return x * tan(x) - 1;
}

// Obliczanie kwadratur złożonych dla funkcji jednej zmiennej

double quad_rect(Func1vFp f1, double a, double b, double c)
{
	return (b - a) * f1(c);
}

double lerp(double a, double b, double r)
{
	return b * r + a * (1 - r);
}

double quad_rect_left(Func1vFp f1, double a, double b, int n) {  // Prostokatow leftpoint 
	double s = 0;
	double sp = a; // start point
	double ep = a; // end point
	
	for (int i = 1; i <= n; i++)
	{
		sp = ep;
		ep = lerp(a, b, (double)i / n);
		s += quad_rect(f1, sp, ep, sp);	
	}
	return s;
}

double quad_rect_right(Func1vFp f1, double a, double b, int n) {  // Prostokatow rightpoint 
	double s = 0;
	double sp = a; // start point
	double ep = a; // end point
	
	for (int i = 1; i <= n; i++)
	{
		sp = ep;
		ep = lerp(a, b, (double)i / n);
		s += quad_rect(f1, sp, ep, ep);	
	}
	return s;
}

double quad_rect_mid(Func1vFp f1, double a, double b, int n) {  // Prostokatow midpoint 
	double s = 0;
	double sp = a; // start point
	double ep = a; // end point
	
	for (int i = 1; i <= n; i++)
	{
		sp = ep;
		ep = lerp(a, b, (double)i / n);
		s += quad_rect(f1, sp, ep, (sp + ep) / 2);	
	}
	return s;
}


double quad_trap(Func1vFp f1, double a, double b, int n) {  // Trapezow
	double s = 0;
	double sp = a; // start point
	double ep = a; // end point
	
	for (int i = 1; i <= n; i++)
	{
		sp = ep;
		ep = lerp(a, b, (double)i / n);
		s += (ep - sp) / 2 * (f1(sp) + f1(ep));	
	}
	return s;
}

double quad_simpson(Func1vFp f1, double a, double b, int n) {  // Simpsona
	double s = 0;
	double sp = a; // start point
	double ep = a; // end point
	
	for (int i = 1; i <= n; i++)
	{
		sp = ep;
		ep = lerp(a, b, (double)i / n);
		s += (ep - sp) / 6 * (f1(sp) + 4 * f1((sp + ep) / 2) + f1(ep));	
	}
	return s;
}

// Definiowanie nazwy dla typu wskaznika do funkcji obliczającej kwadraturę funkcji jednej zmiennej  
typedef double (*QuadratureFp)(Func1vFp, double, double, int);  

// Definicja i inicjowanie tablicy wskaznikow do funkcji jednej zmiennej, w kolejności: f_poly, f_rat, f_trig, f_exp.
Func1vFp func_tab[4]={f_poly, f_rat, f_trig, f_exp};
// Definicja i inicjowanie tablicy wskaznikow do kwadratur dla funkcji jednej zmiennej, w kolejności:
// quad_rect_left, quad_rect_right, quad_rect_mid, quad_trap, quad_simpson.
QuadratureFp quad_tab[5]={quad_rect_left, quad_rect_right, quad_rect_mid, quad_trap, quad_simpson};

// Oblicza i zwraca wartosc wskazanej indeksem quad_no kwadratury dla wskazanej indeksem fun_no funkcji
// w przedziale [a,b] z podzialem na n podprzedzialow. 
double quad_select(int func_no, int quad_no, double a, double b, int n) {
	return quad_tab[quad_no](func_tab[func_no], a, b, n);
}

// Algorytm adaptacyjny obliczania kwadratury,
double recurs(Func1vFp f, double a, double b, double S, double delta, QuadratureFp quad, int level) {
	double ss = quad(f, a, b, 2);
	if (ss <= S + delta && ss >= S - delta) return ss;
	if (level > RECURS_LEVEL_MAX) return nan("");
	double c = (a + b) / 2;
	return recurs(f, a, c, S, delta / 2, quad, level + 1) + 
		recurs(f, c, b, S, delta / 2, quad, level + 1);
}

// Funkcja inicjująca rekurencję
double init_recurs(Func1vFp f, double a, double b, double delta, QuadratureFp quad) {
	return recurs(f, a, b, 0, delta, quad, 0); 
}

///////////////////////////////////////////////////////////////
// 7.2 Calki dwukrotne:

typedef double (*Func2vFp)(double, double);  // Definiowanie nazwy dla typu wskaznika do funkcji dwoch zmiennych

// Przykladowe funkcje dwoch zmiennych:
double func2v_1(double x, double y) {  //fxy1
    return sin(x/y);
}
double func2v_2(double x, double y) {
    return 2-x*x-y*y*y;
}

// Przykladowe funkcje brzegu obszaru calkowania (do calki podwojnej nad obszarem normalnym) 
double lower_bound1(double x) {
    return sin(x);
}
double upper_bound1(double x) {
    return x*x+1;
}
double lower_bound2(double x) {
    return 0.7*exp(-2*x*x);
}
double upper_bound2(double x) {
    return sin(10*x);
}

// Metoda prostokatow (leftpoint) oblicza przyblizenie calki podwojnej nad obszarem prostokatnym
double dbl_integr(Func2vFp f, double x1, double x2, int nx, double y1, double y2, int ny)  {
}

// Oblicza kwadrature prostokatow midpoint dla calki podwojnej nad obszarem normalnym wzgledem osi 0x
double dbl_integr_normal_1(Func2vFp f, double x1, double x2, int nx, double hy, Func1vFp fg, Func1vFp fh)  {
}

// Oblicza kwadrature prostokatow leftpoint dla calki podwojnej nad obszarami normalnymi wzgledem osi 0x 
double dbl_integr_normal_n(Func2vFp f, double x1, double x2, int nx, double y1, double y2, int ny, Func1vFp fg, Func1vFp fh)  {
}

///////////////////////////////////////////////////////////////
// 7.3 Calki wielokrotne:

typedef double (*FuncNvFp)(const double*, int);
typedef int (*BoundNvFp)(const double*, int);

// Przykladowa funkcja trzech zmiennych:
double func3v(const double v[], int n) {
    return v[0]-v[1]+2*v[2];
}

// Przykładowy predykat w przestrzeni 3-wymiarowej
int bound3v(const double v[], int n) {
    if(v[0]>0 && v[0]<0.5 && v[1]*v[1]+(v[2]-1)*(v[2]-1)<1) return 1;   // walec
    return 0;
}

// Przykladowa funkcja n zmiennych:
double funcNv(const double v[], int n) {
    double fv=1.;
    for(int i=1; i<n; ++i) fv += sin(i*v[i]);   
    return fv;
}
// Przykładowy predykat w przestrzeni n-wymiarowej
int boundNv(const double v[], int n) {
    double r=0.0;
    for(int i=0; i<n; ++i) r += (v[i]-1)*(v[i]-1);  // hiperkula n-wymiarowa
    if(r > 1.) return 0;
    return 1;
}

// Obliczanie calek wielokrotnych

// Oblicza calke potrojna "nad" prostopadloscianem z predykatem wykluczajacym jego czesci (jezeli boundary != NULL).
// Metoda prostokatow wstecz (rightpoint) wzdluz kazdej zmiennej.
double trpl_quad_rect(FuncNvFp f, const double variable_lim[][2], const int tn[], BoundNvFp boundary) {     
}

// Oblicza calke wielokrotna (funkcji n zmiennych) "nad" n wymiarowym hiperprostopadloscianem z predykatem wykluczajacym jego czesci (jezeli boundary != NULL).
// Metoda prostokatow midpoint wzdluz kazdej zmiennej.
void recur_quad_rect_mid(double *psum, FuncNvFp f, int variable_no, double tvariable[], const double variable_lim[][2], const int tn[], int level, BoundNvFp boundary) {
}

int main(void)
{
    int to_do, n, nx, ny, integrand_fun_no, method_no, flag;
    double a,b,x1,x2,y1,y2,hy,sum,delta;
    double tvariable[N_MAX], variable_lim[N_MAX][2];
    int tn[N_MAX];

    if (TEST) printf("Wpisz numer testu [1, 7]: ");
    scanf("%d", &to_do);
    switch (to_do) {
    case 1: // 7.1.1 wybor funkcji i metody
        if (TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow: ");
        scanf("%lf %lf %d", &a, &b, &n);
        for(int q=0; q<5; ++q) {
            for(int f=0; f<4; ++f) {
                    printf("%.5f ",quad_select(f, q, a, b, n));
            }
            printf("\n");
        }
        break;
    case 2: // 7.1.2 rekurencyjny algorytm adaptacyjny
        if(TEST) printf("Nr funkcji (0-3) i metody (0-4): ");
        scanf("%d %d",&integrand_fun_no,&method_no);
        if (TEST) printf("Wpisz przedzial calkowania i tolerancje bledu: ");
        scanf("%lf %lf %lf", &a, &b, &delta);
        printf("%.5f\n",init_recurs(func_tab[integrand_fun_no],a,b,delta,quad_tab[method_no]));
        break;
    case 3: // 7.2.1 calka podwojna nad prostokatem
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz y: ");
        scanf("%lf %lf %d",&y1,&y2,&ny);
        printf("%.5f\n",dbl_integr(func2v_2, x1, x2, nx, y1, y2, ny));
        break;
    case 4: // 7.2.2 calka podwojna nad obszarem normalnym
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow zmiennej x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz dlugosc podprzedzialu wzdluz y: ");
        scanf("%lf",&hy);
        printf("%.5f\n",dbl_integr_normal_1(func2v_2, x1, x2, nx, hy, lower_bound2, upper_bound2));
        break;
    case 5: // 7.2.3 calka podwojna nad wieloma obszarami normalnymi
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz y: ");
        scanf("%lf %lf %d",&y1,&y2,&ny);
        printf("%.5f\n",dbl_integr_normal_n(func2v_2, x1, x2, nx, y1, y2, ny, lower_bound2, upper_bound2));
        break;
    case 6: // 7.3.1 calka potrojna po prostopadloscianie
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 1. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[0][0],&variable_lim[0][1],tn);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 2. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[1][0],&variable_lim[1][1],tn+1);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 3. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[2][0],&variable_lim[2][1],tn+2);
        if(TEST) printf("Wpisz 1 gdy ograniczenie obszaru całkowania ma byc aktywne; 0 - w przeciwnym przypadku: ");
        scanf("%d",&flag);
        printf("%.5f\n",trpl_quad_rect(func3v, variable_lim, tn, flag?bound3v:NULL));
        break;
    case 7: // 7.3.3 calka n-wymiarowa na hiperprostopadloscianie
        if(TEST) printf("Wpisz liczbe zmiennych <= %d: ", N_MAX);
        scanf("%d",&n);
        if(n>N_MAX) break;
        for(int i=0; i<n; ++i) {
            if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow %d. zmiennej: ",i+1);
            scanf("%lf %lf %d",&variable_lim[i][0],&variable_lim[i][1],tn+i);
        }
        if(TEST) printf("Wpisz 1 gdy ograniczenie obszaru całkowania ma byc aktywne; 0 - w przeciwnym przypadku: ");
        scanf("%d",&flag);
        sum=0.;
        recur_quad_rect_mid(&sum, funcNv, n, tvariable, variable_lim, tn, 0, flag?boundNv:NULL);//        recur_quad_rect_mid(&sum, func3v_1, n, tvariable, variable_lim, tn, 0, flag?bound3v_2:NULL);
//        recur_quad_rect(&sum, funcNv, n, tvariable, variable_lim, tn, tqr, 0, flag?boundNv:NUL
        printf("%.5f\n", sum);
	break;
    default:
        printf("Numer testu spoza zakresu [1, 7] %d", to_do);
    }
    return 0;
}


