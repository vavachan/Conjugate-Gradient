extern int N;
extern  double BOX;
extern  double epsilon;
extern  double magnitude;
extern  double p_applied;

extern  double *ATOMx;
extern  double *ATOMy;
extern  double *ATOMz;

extern  double *tATOMx;
extern  double *tATOMy;
//extern  double *tATOMz;

extern  double *RAD  ;
extern  double *CONJ_D;
extern  double *G;
extern  double *X;
extern  double delta ;
double U_ij( double , double );

double der_U_ij( double , double );

double energy( double *,int );

double energy_force( double *,  double *, int);
void calculate_gradient( double *,  double *, int);

void write_config(double ,  double *, char*);

void make_list();
