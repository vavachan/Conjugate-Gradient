 double function1dim( double , double *, double *,  double (*f)( double *, int),int);

void bracket( double &,  double &,  double &,  double &,  double &,  double &,  double *, double *,  double (*f) ( double *,int),  int );

 double brent( double ,  double ,  double ,  double ,  double ,  double , double &,  double *,  double *,  double (*f)( double *, int), int );

void line_search( double &,  double *, double *,  double (*fdf)( double *, double*, int), int );

double minimize( int ,  double *, int ,  double (*df)( double *,  double *, int), void (*ws)(double ,  double *, char*));
