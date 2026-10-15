/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Copyright (C) 2006 G. Ziegenhain

    Nichtergodizitaetsparameter (NEP) der Modenkopplungstheorie fuer
    monodisperse Systeme in zwei Dimensionen (Listing der Diplomarbeit 2006).

    2026: Schalter QUAD
      QUAD=1 (Standard): exakte Zellintegration der integrierbaren
             1/sqrt-Singularitaeten des 2D-Vertex bei p=|q-k| und p=q+k.
             Ergibt eta_c^PY = 0.6915 (Literatur: 0.6914, Hajnal et al. 2009).
      QUAD=0: Originalschema 2006 (Mittelpunkte + Simpson-Randgewichte,
             letzte Zelle fehlt) -> eta_c^PY = 0.7044 der Arbeit.
    Startwert jetzt f(q) = 1 (STARTVALUE 3); STARTVALUE 2 liest 'fq_infile'.

    Aufruf: nep <punkte> <dichte n*sigma^2> <qmax> <sq-datei> <zeilen der sq-datei>
    Ausgabe: Datei 'f' + <sq-datei ohne erstes Zeichen>; Exitcode 0 = Glas (|f| >= 1),
             sonst Fluessigkeit.
*/
#include <stdio.h>
#include <stdlib.h>
#define __USE_GNU
#define __USE_ISOC99
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

//#define __DEBUG__

#ifndef STARTVALUE
#define STARTVALUE 3
#endif
#ifndef QUAD
#define QUAD 1
#endif
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L
#endif



typedef struct {
	double *matrix;	
	double *ptr, *ptr1;
	unsigned int size;
} t_matrix;

t_matrix init_matrix (unsigned int size) {
	t_matrix matrix = {
		.size = size,
		.matrix = (double *) calloc ((size_t)(size), (size_t)(sizeof (double))),
		.ptr = NULL
	};
	if (!matrix.matrix)
		exit (-1);
	return matrix;
}

void cleanup_matrix (t_matrix matrix) {
	free (matrix.matrix);
}



int out (const char *format, ...) {
        va_list arg;
        int done;
        
        va_start (arg, format);
        done = vfprintf (stdout, format, arg);
        fflush (stdout); 
        va_end (arg);   
                        
        return done;            
}                       
                        
void die (const char *format, ...) { 
        va_list arg;            
        int done;       
                
        va_start (arg, format);
        done = vfprintf (stdout, format, arg);
        fflush (stdout);
        va_end (arg);
                
        exit (-1);
}                       

void matrix_out (t_matrix x, t_matrix y, char *filename) {
	FILE *out_file = fopen (filename, "w");
	unsigned int i;
	
	if (!out_file)
		die ("matrix_out");
	
	for (i = 0, x.ptr = x.matrix, y.ptr = y.matrix; i < x.size; i++, x.ptr++, y.ptr++)
		fprintf (out_file, "\n%lf %lf", *x.ptr, *y.ptr);
	
	fflush (out_file);
	fclose (out_file);
}

void read_function (t_matrix x, t_matrix function, char *name) {
	unsigned int i;
	FILE *function_file = fopen (name, "r");
	if (!function_file)
		die ("\nread_function: %s\n", name);

	for (i = 0, x.ptr = x.matrix, function.ptr = function.matrix; i < function.size; i++, x.ptr++, function.ptr++) 
		fscanf (function_file, "\n%lf %lf", x.ptr, function.ptr);
	
	fclose (function_file);
}



void spline (t_matrix x, t_matrix y, t_matrix y2) {
	unsigned int i, k;
	double p, qn, sig, un;
	t_matrix u = init_matrix (x.size-1);

	u.matrix[0] = 0.0;

	for (i = 1; i < x.size-1; i++) {
		sig = (x.matrix[i]-x.matrix[i-1])/(x.matrix[i+1]-x.matrix[i-1]);
		p = sig*y2.matrix[i-1]+2.0;
		y2.matrix[i] = (sig-1.0)/p;
		u.matrix[i] = (y.matrix[i+1]-y.matrix[i])/(x.matrix[i+1]-x.matrix[i])-(y.matrix[i]-y.matrix[i-1])/(x.matrix[i]-x.matrix[i-1]);
		u.matrix[i] = (6.0*u.matrix[i]/(x.matrix[i+1]-x.matrix[i-1])-sig*u.matrix[i-1])/p;
	}

	qn = un = 0.0;

	y2.matrix[y2.size-1] = (un-qn*u.matrix[y2.size-2])/(qn*y2.matrix[y2.size-2]+1.0);
	for (k = y2.size-2; k > 0; k--) 
		y2.matrix[k] = y2.matrix[k]*y2.matrix[k+1]+u.matrix[k];

	cleanup_matrix (u);
}

double splint (t_matrix x, t_matrix y, t_matrix y2, double x_) {
	int klo = 0, khi = x.size-1, k;
	double h, b, a;

	while (khi-klo > 1) {
		k = (khi+klo) >> 1;
		if (x.matrix[k] > x_) 
			khi = k;
		else
			klo = k;
	}

	h = x.matrix[khi]-x.matrix[klo];
	if (h == 0.0)
		return (double)(0.0);
	a = (x.matrix[khi]-x_)/h;
	b = (x_-x.matrix[klo])/h;
	return (double)(a*y.matrix[klo]+b*y.matrix[khi]+((a*a*a-a)*y2.matrix[klo]+(b*b*b-b)*y2.matrix[khi])*(h*h)/6.0);
}



int main (int argc, char *argv[]) {
	out ("\n *** NEP monodispers ***");

	out ("\n-> Init");
	unsigned int l,		
		     i, j;
	unsigned long int itera = 0;

	unsigned int size = 300, in_size = 1020;
	double n = 0.6;
	double qmin = 0.5, qmax = 50.0;
	char *sq_infilename = NULL;
	out ("\nNur N, qmax angeben ( Sperls Diskretisierung )");
	if (argc != 1+5)
		die ("\ncmd <points> <density [s^2]> <qmax [s]> <sq_infilename> <insize>\n");
	else {
		size = (unsigned int)(atoi (argv[1]));
		n = (double)(atof (argv[2]));
		qmax = (double)(atof (argv[3]));
		sq_infilename = argv[4];
		in_size = (unsigned int)(atoi (argv[5]));
		qmin = .5*qmax/(double)(size);
	}

	out ("\n\t* size=%d n=%g qmin=%g qmax=%g sq_infilename:%s in_size:%d", size, n, qmin, qmax, sq_infilename, in_size);
	
	t_matrix sq_bessel = init_matrix (in_size);
	t_matrix cq_bessel = init_matrix (in_size);
	t_matrix q_bessel = init_matrix (in_size);
	t_matrix sq_deriv = init_matrix (in_size);
	t_matrix cq_deriv = init_matrix (in_size);

	t_matrix sq = init_matrix (size);
	t_matrix cq = init_matrix (size);

	t_matrix q = init_matrix (size);
	t_matrix fq_new = init_matrix (size);
	t_matrix fq_old = init_matrix (size);
	t_matrix k = init_matrix (size);
	t_matrix p = init_matrix (size);
	t_matrix mct = init_matrix (size);
	

	out ("\n\t* qkp");
	double dq = qmax/(double)(size);
	for (i = 0, q.ptr = q.matrix, p.ptr = p.matrix, k.ptr = k.matrix; i < q.size; i++, q.ptr++, p.ptr++, k.ptr++)
		*q.ptr = *p.ptr = *k.ptr = qmin+(double)(i)*dq;

	out ("\n\t* S(q)");
	FILE *sq_infile = fopen (sq_infilename, "r");
	if (!sq_infile)	
		die ("sq_infile");
	for (i = 0, sq_bessel.ptr = sq_bessel.matrix, q_bessel.ptr = q_bessel.matrix; i < sq_bessel.size; i++, sq_bessel.ptr++, q_bessel.ptr++) 
		fscanf (sq_infile, "\n%lf %lf", q_bessel.ptr, sq_bessel.ptr);
	fclose (sq_infile);
	if (qmin < q_bessel.matrix[0] || qmax > q_bessel.matrix[q_bessel.size-1])
		out ("\nqmin:%lf %lf qmax: %lf %lf\n", qmin, q_bessel.matrix[0], qmax, q_bessel.matrix[q_bessel.size-1]);
	spline (q_bessel, sq_bessel, sq_deriv);
	for (i = 0, sq.ptr = sq.matrix, q.ptr = q.matrix; i < sq.size; i++, sq.ptr++, q.ptr++) 
		*sq.ptr = splint (q_bessel, sq_bessel, sq_deriv, *q.ptr);
	
	out (" c(q)");
	for (i = 0, cq.ptr = cq.matrix, sq.ptr = sq.matrix; i < cq.size; i++, cq.ptr++, sq.ptr++)
		*cq.ptr = (*sq.ptr-1.0)/ *sq.ptr / n;
	
	
	out ("\n\t* Iterationskonstanten");
	double iconst = n/8.0*(double)(powl(2.0*M_PIl,-2.0)*dq*dq)*2.0;
	out ("\niconst:%lf", iconst);
	
	unsigned int jmax = 0;
	unsigned int jmin = 0;
	int j_ = 0;
	const unsigned int j_const = (unsigned int)(q.matrix[0]/dq);
	
	double q2, k2, p2;
	double q5;
	double abs_sq;	
	double sq_rt;	
	const double renorm1 = 3.0/8.0, renorm2 = 7.0/6.0, renorm3 = 23.0/24.0;	
	double sqa;	
	double summand_inner = 0.0, summand_outer = 0.0;	
	double step = 0.5;	
	const double step_multiplicator = 1.1;	
	
	double *sp, *sk;
	
	out ("\n-> Startwert f(q) ");
#if STARTVALUE == 1
	out ("f(q) = S(q)");
	for (i = 0, fq_new.ptr = fq_new.matrix, sq.ptr = sq.matrix; i < fq_new.size; i++, fq_new.ptr++, sq.ptr++)
		*fq_new.ptr = *sq.ptr;
#elif STARTVALUE == 3
	out ("f(q) = 1");
	for (i = 0; i < fq_new.size; i++) fq_new.matrix[i] = 1.0;
#elif STARTVALUE == 2	
	out ("einlesen (muß gleiche Skalierung der qs haben!)  'fq_infile'");

	t_matrix q_tmp = init_matrix (q.size);
	read_function (q_tmp, fq_new, "fq_infile"); 
	cleanup_matrix (q_tmp);
	step = 0.9;
	out ("   passe auch Schrittweite an: %lf", step);
#endif
	
	out ("\n-> Beginne Iteration");
	double norm_fq_new = 0.0, norm_fq_old = 0.0, norm_fq_adap1 = 0.0;
	double norm_fq_max = 1e-8;
	
	time_t start, now;	
        time (&start);

	do {
		norm_fq_adap1 = norm_fq_old;
		norm_fq_old = norm_fq_new;
		norm_fq_new = 0.0;
		for (i = 0, fq_new.ptr = fq_new.matrix, fq_old.ptr = fq_old.matrix; i < fq_new.size; i++, fq_new.ptr++, fq_old.ptr++) {
			*fq_old.ptr = *fq_new.ptr;
			norm_fq_new += *fq_new.ptr**fq_new.ptr;
		}
		norm_fq_new = sqrt (norm_fq_new);
		out ("\n(i:%ld) (n:%.8f) (d:%.8f) (s:%f)", ++itera, norm_fq_new, fabs(norm_fq_old-norm_fq_new), step);

		out (" m");
		for (l = 0, mct.ptr = mct.matrix, q.ptr = q.matrix, sq.ptr = sq.matrix; l < mct.size; l++, mct.ptr++, q.ptr++, sq.ptr++) {
			*mct.ptr = 0.0;
			q5 = 1.0/(*q.ptr**q.ptr**q.ptr**q.ptr**q.ptr);
			for (i = 0, k.ptr = k.matrix, sk = sq.matrix, fq_new.ptr = fq_new.matrix, cq.ptr = cq.matrix; i < sq.size; i++, k.ptr++, sk++, fq_new.ptr++, cq.ptr++) {
				summand_outer = 0.0;
				j_ = abs(l-i)-(int)(j_const);
				jmin = (j_ > 0) ? (unsigned int)(j_) : 0;
#if QUAD == 1
				j_ = i+l+1;
#else
				j_ = i+l+(int)(j_const);
#endif
				jmax = ((unsigned int)(j_) < size) ? j_ : (unsigned int)(size);
				for (j = jmin, p.ptr = p.matrix+jmin, sp = sq.matrix+jmin, fq_new.ptr1 = fq_new.matrix+jmin, cq.ptr1 = cq.matrix+jmin; j < jmax; j++, p.ptr++, sp++, fq_new.ptr1++, cq.ptr1++) {
					q2 = *q.ptr**q.ptr;
					p2 = *p.ptr**p.ptr;
					k2 = *k.ptr**k.ptr;
					abs_sq = (q2-p2+k2)**cq.ptr+(q2-k2+p2)**cq.ptr1;
					abs_sq *= abs_sq;
#if QUAD == 1
					{	/* exakte Zellintegration: int_a^b p dp / sqrt(k^2-((q^2+k^2-p^2)/2q)^2) = q [asin(w_a/A)-asin(w_b/A)] */
						double A = 2.0**q.ptr**k.ptr, pa = *p.ptr-0.5*dq, pb = *p.ptr+0.5*dq;
						double wa = (q2+k2-pa*pa)/A, wb = (q2+k2-pb*pb)/A;
						wa = (wa > 1.0) ? 1.0 : ((wa < -1.0) ? -1.0 : wa);
						wb = (wb > 1.0) ? 1.0 : ((wb < -1.0) ? -1.0 : wb);
						sqa = *q.ptr*(asin(wa)-asin(wb))/(dq**p.ptr);
					}
#else
					sq_rt = (q2-p2+k2)/(2.0**q.ptr);
					sq_rt *= sq_rt;
					sqa = pow(k2-sq_rt,-0.5);
#endif

					
					if (isnan (sqa) || isinf (sqa)) 
						die ("nan/inf: inner");
					
					summand_inner = *p.ptr**sp*sqa*abs_sq**fq_new.ptr1;

#if QUAD == 0
					if (j == 0 || j == jmax-1) 
						summand_inner *= renorm1;
					else if (j == 1 || j == jmax-2) 
						summand_inner *= renorm2;
					else if (j == 2 || j == jmax-3) 
						summand_inner *= renorm3;
#endif
						
					summand_outer += summand_inner;
				}
				summand_outer *= *k.ptr**sk**fq_new.ptr*q5;

				if (isnan (summand_outer) || isinf (summand_outer)) 
					die ("nan/inf: outer");

#if QUAD == 0
				if (i == 0 || i == sq.size-1) 
					summand_outer *= renorm1;
				else if (i == 1 || i == sq.size-2) 
					summand_outer *= renorm2;
				else if (i == 2 || i == sq.size-3) 
					summand_outer *= renorm3;
#endif
		
				*mct.ptr += summand_outer;
			}
			*mct.ptr *= iconst**sq.ptr;
			
			if (isnan (*mct.ptr) || isinf (*mct.ptr)) 
				die ("nan/inf: mct");
		}

		if (fabs(norm_fq_adap1-norm_fq_old) > fabs(norm_fq_old-norm_fq_new) && step*step_multiplicator < 0.999)
			step *= step_multiplicator;
		
		out ("n");
		for (i = 0, fq_new.ptr = fq_new.matrix, fq_old.ptr = fq_old.matrix, mct.ptr = mct.matrix; i < fq_new.size; i++, fq_new.ptr++, fq_old.ptr++, mct.ptr++)
			*fq_new.ptr = step*(*mct.ptr/(1.0+*mct.ptr)) + (1.0-step)**fq_old.ptr;
		
		time (&now);
		out (" time: %lf", difftime (now, start));
	} while (fabs(norm_fq_old-norm_fq_new) > norm_fq_max || step < 0.9);

	
	goto out;
out:
	out ("\n-> f(q) ausgeben");
        char *tmp = (char *) &(*(sq_infilename+1));
	char fq_outfilename[255];
	sprintf (fq_outfilename, "f%s", tmp);
			
	FILE *fq_outfile = fopen (fq_outfilename, "w");
 	for (i = 0, fq_new.ptr = fq_new.matrix, q.ptr = q.matrix; i < fq_new.size; i++, fq_new.ptr++, q.ptr++)
		fprintf (fq_outfile, "\n%.32lf %.32lf", *q.ptr, *fq_new.ptr);
	fclose (fq_outfile);

	
	goto clean;	
clean:
	out ("\n-> Cleanup");
	
	cleanup_matrix (sq);
	cleanup_matrix (cq);

	cleanup_matrix (sq_bessel);
	cleanup_matrix (cq_bessel);
	cleanup_matrix (sq_deriv);
	cleanup_matrix (cq_deriv);
	cleanup_matrix (q_bessel);

	cleanup_matrix (q);
	cleanup_matrix (fq_new);
	cleanup_matrix (fq_old);
	cleanup_matrix (k);
	cleanup_matrix (p);
	cleanup_matrix (mct);
	
	out ("\nok.");
	if (norm_fq_new < 1.0) {
		exit (-1);
	} 
	return (0);
}
