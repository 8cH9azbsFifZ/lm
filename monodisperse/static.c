/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2007 G. Ziegenhain, gerolf@ziegenhain.com
*/

#include <stdio.h>
#include <stdlib.h>
#define __USE_GNU 1 // für M_PIl: long double PI aus libmath
#define _GNU_SOURCE // Für j0l: long double Bessel aus libmath
#include <math.h>
long double j0l(long double x);
long double j1l(long double x);
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <time.h>


#define __USE_GSL 1 // GNU Scietific Library - sonst genäherte Besselnullstellen
#ifdef __USE_GSL
#include <gsl/gsl_sf_bessel.h>
#else   // Genäherte Formel, die ab 2. Nullstelle auf mind. 5 Stellen genau ist
long double gsl_sf_bessel_zero_J0 (int number) {
	return M_PIl/4.0*(4.0*number-1.0) +
		1.0/(2.0*M_PIl*(4.0*number-1.0)) -
		31.0/(6.0*powl(M_PIl*(4.0*number-1.0),3.0)) +
		3779.0/(15.0*powl(M_PIl*(4.0*number-1.0),5.0));
}
#endif

#define NDEBUG   // Für assert
//#define DEBUG	// Funktionenausgabe


/******************************
 * Schalter für den Algorithmus
 ******************************/
#define SYSTEM 2
/* 1: Harte Scheibchen
 * 2: repulsiv dipolare harte Scheibchen
 * 3: LJ
 * 4: attraktiv dipolare harte Scheibchen
 */

// Dichte: einzugeben in n*sigma^2
// phi = n*sigma^2 * pi/4 in 2d
// für vgl in 3d: phi_3d = n*sigma^2 * pi/6

#define APPROXIMATION 2   	// Zu verwendende Approximation:
				// 1: PY 
				// 2: HNC
#define STARTVALUE 1   	// Startwert:
			// 1: Gamma(r) == 0
			// 2: Gamma(r) == "Fitfunktion"
			// 3: Gamma(r) aus "startvalue" einlesen

//#define ASYMPTOTE 2
/* Asymptoten
 * 1: r < 1: H(r) = -1
 * 2: r < pow(n,-.5): H(r) = -1
 */

//#define DYN_COARSE
//#define DYN_COARSE_2
//#define PIC_AFTER
//#define DYN_BASIS
//#define DYN_BASIS_1
//#define NEW_BASIS
//#define SLOW_PICARD
#define BASIS_LINE

/* Testroutinen:
 */
//#define TEST 1


/******************* 
 * Ausgabefunktionen
 *******************/
/* Zeitechte Konsolenausgabe */
int out (const char *format, ...) {
        va_list arg;
        int done;
        
        va_start (arg, format);
        done = vfprintf (stdout, format, arg);
        fflush (stdout); 
        va_end (arg);   
                        
        return done;            
}                       
                        
/* Programm mit Fehlermeldung beenden */
void die (const char *format, ...) { 
        va_list arg;            
        int done;       
                
        va_start (arg, format);
        done = vfprintf (stdout, format, arg);
        fflush (stdout);
        va_end (arg);
                
        exit (-1);
}                       


/*******************************************
 * Allokieren von Speicherplatz für Matrizen
 *******************************************/

typedef struct {
	long double *matrix;		// Start des Arrays
	long double *ptr, *ptr1;	// Zeiger in dem Array
	unsigned int size;	// Größe des Arrays
} t_matrix;

t_matrix init_matrix (unsigned int size) {
	t_matrix matrix = {
		.size = size,
		.matrix = (long double *) calloc ((size_t)(size), (size_t)(sizeof (long double))),
		.ptr = NULL
	};
	if (!matrix.matrix)
		exit (-1);
	return matrix;
}

void cleanup_matrix (t_matrix matrix) {
	free (matrix.matrix);
}


long double *alloc_matrix_1 (int anz_1) {   // Liefert einen Vektor mit Einträgen = 0
    	long double *m;
    	int i;
    
    	m = (long double *) malloc ((size_t) (anz_1 * sizeof (long double)));
    	if (!m) 
		die ("allocation failure in alloc_matrix_1");

    	for (i = 0; i < anz_1; i++) 
	    	m[i] = 0.0;
    
    	return m;
}

long double **alloc_matrix_2(int anz_1,int anz_2) {
    	int i;
    	long double **m;
    
    	m = (long double **) malloc((size_t) (anz_1 * sizeof(long double *)));
    	if (!m) 
		die ("allocation failure in alloc_matrix_2");
    	for (i = 0; i < anz_1; i++) 
		m[i] = alloc_matrix_1(anz_2);
   	return m;
}

void free_matrix_1(long double *m) { 
    	free(m);
}

void free_matrix_2(long double **m,int anz_1) { 
    	int i; 
    
    	for (i = 0; i < anz_1; i++)
		free_matrix_1(m[i]);
    	free(m);
}



/********************
 * Typendeklarationen
 ********************/
/* Struktur für Zerlegung / Komposition */
typedef struct {
	long double *function;
	long double *coarse, *fine;
	long double **P_, **Q_;
	int gitter, nu;
	long double sum;   // temporäre Zwischensumme
} t_composite;

/* Struktur für Konstanten und Fkt. der FBT */
typedef struct {
	long double *F_r, *F_q;
	int gitter;
	long double **bessel;
	long double cons;
} t_fbt;

/* Struktur mit Daten zur Berechnung eines neuen a_m in Marquardtschleife */
typedef struct {
	long double *Cr, *Cq, *Gammar, *Gammaq;
	long double *r_;
	long double p;
	int gitter, nu;
	long double *mayer, *Vr;
	t_fbt fbt1, fbt2;
	t_composite comp, decomp;
} lm_data_type;

/* Debugausgabe von Grob-/Fein-Anteil */
typedef struct {
	long double *coarse;
	long double *fine;
	int coarse_dim, fine_dim;
	int *pic, *nr;
	long double *r_;
	long double **P_;
} t_coarse_fine;

/* Struktur für Ausgabe aller Korrelatoren */
typedef struct {
	long double *Gamma_r, *Gamma_q;
	long double *C_r, *C_q;
	long double *r_, *q_;
	int *gitter;
	int *iteration;
	long double *p;
} t_functions;





/* Funktionen für Marquardt-Algorithmus */
#include "marquardt.c"



/*
 * Interpolation mit Splines
 */
void spline (t_matrix x, t_matrix y, t_matrix y2) {
#define __SPLINE__
	//double yp1 = 0.0, ypn = 0.0;	// 1st derivative at boundary points =0 => natural spline
	unsigned int i, k;
	long double p, qn, sig, un;
#ifdef __SPLINE__
	printf ("=> Init"); fflush (stdout);
	t_matrix u = init_matrix (x.size-1);
#endif

	u.matrix[0] = 0.0;	// Natural spline

#ifdef __SPLINE__
	printf (" for"); fflush (stdout);
#endif
	for (i = 1; i < x.size-1; i++) {
		sig = (x.matrix[i]-x.matrix[i-1])/(x.matrix[i+1]-x.matrix[i-1]);
		p = sig*y2.matrix[i-1]+2.0;
		y2.matrix[i] = (sig-1.0)/p;
		u.matrix[i] = (y.matrix[i+1]-y.matrix[i])/(x.matrix[i+1]-x.matrix[i])-(y.matrix[i]-y.matrix[i-1])/(x.matrix[i]-x.matrix[i-1]);
		u.matrix[i] = (6.0*u.matrix[i]/(x.matrix[i+1]-x.matrix[i-1])-sig*u.matrix[i-1])/p;
	}

	qn = un = 0.0;	// Natural spline

#ifdef __SPLINE__
	printf (" for"); fflush (stdout);
#endif
	y2.matrix[y2.size-1] = (un-qn*u.matrix[y2.size-2])/(qn*y2.matrix[y2.size-2]+1.0);
	for (k = y2.size-2; k > 0; k--) 
		y2.matrix[k] = y2.matrix[k]*y2.matrix[k+1]+u.matrix[k];

#ifdef __SPLINE__
	printf (" Cleanup"); fflush (stdout);
#endif
	cleanup_matrix (u);
}

long double splint (t_matrix x, t_matrix y, t_matrix y2, long double x_) {
	int klo = 0, khi = x.size-1, k;
	long double h, b, a;

	while (khi-klo > 1) {
		k = (khi+klo) >> 1;
		if (x.matrix[k] > x_) 
			khi = k;
		else
			klo = k;
	}

	h = x.matrix[khi]-x.matrix[klo];
	if (h == 0.0)
		return (long double)(0.0);
	a = (x.matrix[khi]-x_)/h;
	b = (x_-x.matrix[klo])/h;
	return (long double)(a*y.matrix[klo]+b*y.matrix[khi]+((a*a*a-a)*y2.matrix[klo]+(b*b*b-b)*y2.matrix[khi])*(h*h)/6.0);
}




/***************************************
 * Routinen für die Picard-Hauptschleife
 ***************************************/
/* PY - Approximation */
void PY (long double *C_r, long double *Gamma_r, long double *mayer, int gitter) {
	unsigned int i;

	for (i = 0; i < gitter; i++) 
		C_r[i] = mayer[i] * (1.0+Gamma_r[i]);
}

/* HNC - Approximation */
void HNC (long double *C_r, long double *Gamma_r, long double *Vr, int gitter) {
	unsigned int i;

	for (i = 0; i < gitter; i++) 
		C_r[i] = expl (Vr[i] + Gamma_r[i])-Gamma_r[i]-1.0;
}

/* FT: Orts -> Impulsraum */
void FBT1 (t_fbt fbt1) {
	unsigned int m, i;
	
	for (m = 0; m < fbt1.gitter; m++) {
		fbt1.F_q[m] = 0.0;
		for (i = 0; i < fbt1.gitter; i++) 
			fbt1.F_q[m] += (long double) (fbt1.F_r[i] * fbt1.bessel[m][i]);
		fbt1.F_q[m] *= fbt1.cons;
	}
}

/* OZ - Gleichung */
void OZ (long double *C_q, long double *Gamma_q, long double p, int gitter) {
	unsigned int i;

	for (i = 0; i < gitter; i++) 
		Gamma_q[i] = C_q[i]*p*C_q[i]/(1.0-p*C_q[i]);
}

/* FT: Impulsraum -> Ortsraum */
void FBT2 (t_fbt fbt2) {
	unsigned int i, m;

	for (i = 0; i < fbt2.gitter; i++) {
		fbt2.F_r[i] = 0.0;
		for (m = 0; m < fbt2.gitter; m++) 
			fbt2.F_r[i] += (long double) (fbt2.F_q[m] * fbt2.bessel[m][i]);
		fbt2.F_r[i] *= fbt2.cons;
	}
}



/*******************************************
 * Routinen für die Zerlegungsbasiserzeugung
 *******************************************/
/* Gegeben Pkt. X und Vector r_i liefert dies i: |X-r_i|=minimal */
int find_i (long double *r_, long double X, int gitter) {
	int i = gitter-1, j;
	long double minimum = r_[i];
	
	for (j = 0; j < gitter; j++)
		if (fabsl(r_[j]-X) <= minimum) {
			minimum = fabsl(r_[j]-X);
			i = j;
		}

	return i;
}

/* Invertierung einer 2x2-Matrix */
void invert (long double **matrix, long double **inverse, int dimension) {   
	int k, i, j;   // Counter
	double min_value = powl (10.0, -10.0);

	// Invertierbar hiermit?
	for (i = 0; i < dimension; i++) 
		if (matrix[i][i] < min_value) {
			die ("Invertierfehler: (Diagonalelemet %d) = %.64Lf", i, matrix[i][i]);
			return;
		}

	// Einheitsmatrix
	long double **unity = alloc_matrix_2 (dimension, dimension);
	for (i = 0; i < dimension; i++) for (j = 0; j < dimension; j++) 
		unity[i][j] = (long double)(!(i-j));

	// Invertierung
	long double xmult;
	for (i = 0; i < dimension-1; i++) {
		for (j = i+1; j < dimension; j++) {
			xmult = matrix[j][i]/matrix[i][i];
			for (k = i+1; k < dimension; k++) {
				matrix[j][k] -= xmult*matrix[i][k];
			}
			matrix[j][i] = xmult;
			for (k = 0; k < dimension; k++) {
				unity[j][k] -= xmult*unity[i][k];
			}
		}
	}
	long double sum;
	for (i = 0; i < dimension; i++) {
		inverse[dimension-1][i] = unity[dimension-1][i]/matrix[dimension-1][dimension-1];
		//oder dim-2?
		for (j = dimension-2; j >= 0; j--) {   
			sum = unity[j][i];
			for (k = j+1; k < dimension; k++) {
				sum -= matrix[j][k]*inverse[k][i];
			}
			inverse[j][i] = sum/matrix[j][j];
		}
	}

	free_matrix_2 (unity, dimension);
}

/* Zerlegung in Grob-/Feinanteil */
void decomposite (t_composite decomp) {
//#define _DECOMPOSITE_
	unsigned int i, n;

#ifdef _DECOMPOSITE_
	out ("\n c: n");
#endif
	for (n = 0; n < decomp.nu; n++) 
		decomp.coarse[n] = 0.0;

#ifdef _DECOMPOSITE_
	out (" a");
#endif
	for (n = 0; n < decomp.nu; n++) for (i = 0; i < decomp.gitter; i++) 
		decomp.coarse[n] += decomp.Q_[n][i]*decomp.function[i];

#ifdef _DECOMPOSITE_
	out (" f");
#endif
	for (i = 0; i < decomp.gitter; i++) {
		decomp.sum = 0.0;
		for (n = 0; n < decomp.nu; n++)
			decomp.sum += decomp.coarse[n]*decomp.P_[n][i];
		decomp.fine[i] = decomp.function[i]-decomp.sum;
	}
#ifdef _DECOMPOSITE_
	out (" .");
#endif
}

/* Komposition aus Grob-/Fein-Anteil */
void composite (t_composite comp) {
//#define _COMPOSITE_
	unsigned int i, n;

	for (i = 0; i < comp.gitter; i++) {
#ifdef _COMPOSITE_
		out ("\ni:%d");
#endif
		comp.function[i] = comp.fine[i];
#ifdef _COMPOSITE_
		out (" %Lf", comp.function[i]);
#endif
		for (n = 0; n < comp.nu; n++) {
			comp.function[i] += comp.P_[n][i]*comp.coarse[n];
#ifdef _COMPOSITE_
			out (" n:%d %Lf", n, comp.function[i]);
#endif
		}
	}
}

/* Erzeugung der Basen */
void pq_init (long double *r_i_a, unsigned int nu, long double *r_, unsigned int gitter, long double **P_, long double **Q_) {
//#define _PQ_INIT_	// Frage: warum geht es nicht, wenn die Fkt. ausgegebn werden=????
	unsigned int i, m, n;

	//out ("\n-> Zerlegungsbasis ");
	// Stützpunkte wählen
	int i_[nu];//*i_ = (int *) malloc ((size_t)(nu*sizeof(int)));   // Die nu i_alphas aus den Maximalstellen bestimmen
	for (n = 0; n < nu; n++) { 
		i_[n] = find_i (r_, r_i_a[n], (int)gitter);
#ifdef _PQ_INIT_
		out ("i: %d %Lf", i_[n], r_i_a[n]);
#endif
	}
	
	// *** P_n^i konstruieren ***
	//out ("P_n^i ");
	//*P_ = (long double **) alloc_matrix_2 ((int)nu, (int)gitter);
	for (n = 0; n < nu; n++) 
		if (n == 0) {   // 1. Basis
			for (i = 0; i < gitter; i++) {
				P_[0][i] = 0.0;
				if (i <= i_[0])	P_[0][i] = (long double) (i_[0]-i) / (i_[0]);
			}
		} else if (n == 1) {	// 2. Basis
			for (i = 0; i < gitter; i++) {
				P_[n][i] = 0.0;
				if (0 <= i && i <= i_[1]) P_[n][i] = (long double)(i)/(long double)(i_[1]);
				else if (i_[1] <= i && i <= i_[2]) P_[n][i] = (long double)(i_[2]-i)/(long double)(i_[2]-i_[1]);
			}
		} else {   // Alle anderen Basen
			for (i = 0; i < gitter; i++) {
				P_[n][i] = 0.0;
				if (i_[n-1] <= i && i <= i_[n]) P_[n][i] = fabsl ((long double) (i-i_[n-1]) / (i_[n-1]-i_[n]));
				else if (i_[n] <= i && i <= i_[n+1]) P_[n][i] = fabsl ((long double) (i_[n+1]-i) / (i_[n+1]-i_[n]));
			}
		}

#ifdef _PQ_INIT_
	FILE *p_outfile = fopen ("p_outfile", "w");
	for (n = 0; n < nu; n++)
		for (i = 0; i < gitter; i++)
			fprintf (p_outfile, "\n%d %d %Lf", n, i, P_[n][i]);
	fclose (p_outfile);
#endif
	
	// *** Q_n^i konstruieren ***
	//out ("Q_n^i ");
	//*Q_ = (long double **) alloc_matrix_2 ((int)nu, (int)gitter);
	
	// Dazu die Matrix (PP)_mn = sum_i(P_m^i*P_n^i) berechnen
	long double sum;
	long double **PP = alloc_matrix_2 (nu, nu);
	for (m = 0; m < nu; m++) for (n = 0; n < nu; n++) {
		sum = 0.0;
		for (i = 0; i < gitter; i++)
			sum += P_[m][i]*P_[n][i];
		PP[m][n] = sum;
	}

	// Dazu die Inserve RR von PP berechnen
	long double **RR = alloc_matrix_2 (nu, nu);
	invert (PP, RR, nu);

	// Und schließlich die Qs berechnen
	for (m = 0; m < nu; m++) 
		for (i = 0; i < gitter; i ++) {
			sum = 0.0;
			for (n = 0; n < nu; n++) 
				sum += RR[m][n]*P_[n][i];
			Q_[m][i] = sum;
		}

	//out (" Cleanup");
	//free (i_);
	free_matrix_2 (RR, nu);
	free_matrix_2 (PP, nu);
}



/******************
 * Ausgabe-Routinen
 ******************/
/* Einen Korrelator rausschreiben */
void write_function (long double *function, int gitter, long double *vector, char *name, int iteration) {
	char filename[255];
	sprintf (filename, "%s__%d.dat", name, iteration);
	FILE *function_file = fopen (filename, "w");
	int i;
	
	for (i = 0; i < gitter; i++) 
		fprintf (function_file, "\n%.64Lf   %.64Lf", vector[i], function[i]);
	
	fclose (function_file);
}

/* Grob-/Fein-Anteil rausschreiben */
void coarse_fine_out (t_coarse_fine cf) {
	char filename_coarse[255], filename_fine[255];
	sprintf (filename_coarse, "coarse__pic%d_nr%d", *cf.pic, *cf.nr);
	sprintf (filename_fine, "fine__pic%d_nr%d", *cf.pic, *cf.nr);
	FILE *coarse_file = fopen (filename_coarse, "w"), *fine_file = fopen (filename_fine, "w");

	int i, n;
	
	for (i = 0; i < cf.fine_dim; i++) 
		fprintf (fine_file, "\n%.64Lf %.64Lf", cf.r_[i], cf.fine[i]);
	
	long double tmp;
	for (i = 0; i < cf.fine_dim; i++) {
		fprintf (coarse_file, "\n%.64Lf ", cf.r_[i]);
		tmp = 0.0;
		for (n = 0; n < cf.coarse_dim; n++) 
			tmp += cf.P_[n][i]*cf.coarse[n];
		fprintf (coarse_file, "%.64Lf ", tmp);
	}
	
	fclose (coarse_file);
	fclose (fine_file);
}

/* Rausschreiben aller Korrelatoren */
void function_out (t_functions m) {
	int i;
	long double *tmp = alloc_matrix_1 (*m.gitter);

	// H(r) & H(q)
	for (i = 0; i < *m.gitter; i++)
		tmp[i]= m.Gamma_r[i] + m.C_r[i];
	write_function (tmp, *m.gitter, m.r_, "Hr", *m.iteration);
	for (i = 0; i < *m.gitter; i++)
		tmp[i] = m.Gamma_q[i] + m.C_q[i];
	write_function (tmp, *m.gitter, m.q_, "Hq", *m.iteration);

	// Gamma(r) & Gamma(q)
	write_function (m.Gamma_r, *m.gitter, m.r_, "Gammar", *m.iteration);
	write_function (m.Gamma_q, *m.gitter, m.q_, "Gammaq", *m.iteration);

	// C(r) & C(q)
	write_function (m.C_r, *m.gitter, m.r_, "Cr", *m.iteration);
	write_function (m.C_q, *m.gitter, m.q_, "Cq", *m.iteration);

	// S(q)
	for (i = 0; i < *m.gitter; i++) 
		tmp[i] = 1.0+(*m.p)*(m.Gamma_q[i]+m.C_q[i]);
	write_function (tmp, *m.gitter, m.q_, "Sq", *m.iteration);

	// Aufräumen
	free_matrix_1 (tmp);
}

/* Einen Korrelator einlesen aus Datei */
void read_function (long double *function, long double *r_, int gitter, char *name) {
	unsigned int rows = 0;
	unsigned int i;
	FILE *function_file = fopen (name, "r");
	
	// # Zeilen
	if (!function_file)
		die ("\nfunction_file: %s", name);
	long double tmp;
	do {
		fscanf (function_file, "%Lf %Lf\n", &tmp, &tmp);
		rows++;
	} while (!feof (function_file));
	tmp = 0; // Gegen Warnung
	fclose (function_file);

	// Datein einlesen
	function_file = fopen (name, "r");
	if (!function_file)
		die ("\nfunction_file: %s", name);
	if (gitter < rows) 
		die ("\nKeine Interpolation nötig: %d %d", gitter, rows);
	t_matrix r_in = init_matrix (rows);
	t_matrix f_in = init_matrix (rows);
	t_matrix deriv = init_matrix (rows);
	for (i = 0, r_in.ptr = r_in.matrix, f_in.ptr = f_in.matrix; i < rows; i++, r_in.ptr++, f_in.ptr++) 
		fscanf (function_file, "%Lf %Lf\n", r_in.ptr, f_in.ptr);
	
	// Interpolieren
	spline (r_in, f_in, deriv);
	for (i = 0; i < (unsigned int)gitter; i++) 
		function[i] = splint (r_in, f_in, deriv, r_[i]);

	cleanup_matrix (r_in);
	cleanup_matrix (f_in);
	cleanup_matrix (deriv);
}

/* reduziertes Potential aus experimentellen Parametern zurückgeben */
void return_mayer (long double *mayer, long double *r_, long double *Vr, unsigned int gitter, long double density, long double temperature, long double b_field) {
	unsigned int i;
	
	long double kb = 1.3806505e-23;		// Boltzmannkonstante
	long double chi = 6.2e-12;		// Magnetische Suszeptibilität
	long double mu_ = 1.0e-7;		// mü0/4pi
	long double T = temperature;
	long double B = b_field;
	long double r = 2.35e-6;			// Teilchenradius

	long double sigma = r+r;		// Teilchendurchmesser
	double norm = powl (sigma,-3.0);	// Normierung auf Teilchendurchmesser
	
	long double V_ = -1.0*mu_*chi*chi*B*B/(kb*T)*norm;	// Reduziertes Potential
	out ("\n-> Reduziertes Potential: V:%g", V_);

	// Vr(r) & Mayer(r)
#if SYSTEM == 1	// Harte Scheibchen
	for (i = 0; i < gitter; i++) {
		mayer[i] = (r_[i] < 1.0) ? -1.0 : 0.0;	// Mayerfunktion
		Vr[i] = (r_[i] < 1.0) ? -1.0/0.0 : 0.0;	// reduziertes Potential
	}
#elif SYSTEM == 2	// Dipolare repulsive harte Scheibchen
	for (i = 0; i < gitter; i++) {
		mayer[i] = (r_[i] < 1.0) ? -1.0 : expl (V_*powl(r_[i],-3.0))-1.0;
		Vr[i] = (r_[i] < 1.0) ? -1.0/0.0 : V_*powl(r_[i],-3.0);
	}
#elif SYSTEM == 3	// LJ
	long double epsilon = 1.0;//temperature;

	V_ = -temperature; // ==-1.0*epsilon/(kb*T);
	
	for (i = 0; i < gitter; i++) {
		mayer[i] = (r_[i] < .5) ? -1.0 : expl (V_*(-powl(sigma/r_[i],6.0)+powl(sigma/r_[i],12.0)))-1.0;
		Vr[i] = (r_[i] < .5) ? -1.0/0.0 : V_*(-powl(sigma/r_[i],6.0)+powl(sigma/r_[i],12.0));
	}
#elif SYSTEM == 4	// Dipolare attrative harte Scheibchen
	/* Wichtig: V ~ -2M1M2/r12^3 */
	for (i = 0; i < gitter; i++) {
		mayer[i] = (r_[i] < 1.0) ? -1.0 : expl (-2.0*V_*powl(r_[i],-3.0))-1.0;
		Vr[i] = (r_[i] < 1.0) ? -1.0/0.0 : -2.0*V_*powl(r_[i],-3.0);
	}
#endif
	out (" mittlerer Abstand:%g", pow(density,-.5));
	out (" Vk(r_):%g", Vr[find_i(r_,pow(density,-.5)/kb,(int)gitter)]);
	out ("\n-> Gamma: %g", 10e-7/(kb*T)*B*B*chi*chi*powl(density*powl(sigma,-2.0)*M_PIl,3.0/2.0));
}



/***********
 * Marquardt
 ***********/
/* d=a-a_ zu neuem a berechnen */
void lm_evaluate_default (double* coarse, int nu, double* diff, void *data, int *info) {
//#define _EVAL_
#ifdef _EVAL_
	out ("e");
#endif
	unsigned int i;
	lm_data_type *d = (lm_data_type*)data;

	// Zerlegung
	for (i = 0; i < nu; i++)
		d->comp.coarse[i] = coarse[i];

	composite (d->comp);

	// PY / HNC Approximation
#if APPROXIMATION == 1
	PY (d->Cr, d->Gammar, d->mayer, d->gitter);
#elif APPROXIMATION == 2
	HNC (d->Cr, d->Gammar, d->Vr, d->gitter);
#endif

	// Asymptote für H(r)
#if ASYMPTOTE == 1	// H(r) = -1 r<1
	// <=> gamma(r)+c(r)=h(r)=-1 <=> -1-gamma(r)=c(r)
	for (i = 0; i < d->gitter; i++)
		if (d->r_[i] <= 1.0)
			d->Cr[i] = -1.0-*d->Gammar[i];
		else
			i = d->gitter;
#elif ASYMPTOTE == 2	// H(r) = -1 r<pow(r,-.5)
	for (i = 0; i < d->gitter; i++)
		if (d->r_[i] <= pow (d->p,-.5))
			d->Cr[i] = -1.0-d->Gammar[i];
		else
			i = d->gitter;
#endif

	// OZ-Gleichung
	d->fbt1.F_r = d->Cr; d->fbt1.F_q = d->Cq;
	FBT1 (d->fbt1);
	OZ (d->Cq, d->Gammaq, d->p, d->gitter);
	d->fbt2.F_q = d->Gammaq; d->fbt2.F_r = d->Gammar;
	FBT2 (d->fbt2);
	decomposite (d->decomp);

	for (i = 0; i < nu; i++)
		diff[i] = d->comp.coarse[i]-d->decomp.coarse[i];

#ifdef _EVAL_
	i = 3;
	out ("\n%d c%Lf f%Lf c%Lf f%Lf", i, d->comp.coarse[i], d->comp.fine[i], d->decomp.coarse[i], d->decomp.fine[i]);
#endif
}



/***************
 * Hauptschleife
 ***************/
int main (int argc, char* argv[]) {
	int v;//, w;   // Impuls-Counter
	int i,j;   // Orts-Counter
	int m, n;   // Zerlegungsbasis-Counter / Impuls-Counter (FBT)


	/**********************
	 * Parameter einstellen
	 **********************/
	out ("\n-> Parameter"); fflush (stdout);
	int gitter;				// # Stützpunkte für die Korrelatoren
	long double p;				// reduzierte Dichte: p_=p*sigma^2
	long double Rmax;			// Abschneideradius in 1/sigma^2
	long double temperature, b_field;

#if SYSTEM == 1	// Harte Scheibchen
	out ("\n-> Harte Scheibchen");

#ifndef BASIS_LINE
	if ((int)argc != 1+3)
		die ("\n<cmd> <density> <rmax> <m>\n");
#else
	if ((int)argc <= 1+3)
		die ("\n<cmd> <density> <rmax> <m> [p...]\n");
#endif
	p = (long double) atof (argv[1]);
	gitter = (unsigned int)(atoi(argv[3]));
	Rmax = (unsigned int)(atoi(argv[2]));
	
	temperature = b_field = 0.0;

	out (" p:%Lf Rmax:%Lf gitter:%d", p, Rmax, gitter);
#ifdef BASIS_LINE
	int NU = (int)argc - (1+3);
	int nu = NU;
	long double *r_i_a = (long double*) alloc_matrix_1 ((unsigned int)NU);
	for (i = 0; i < NU; i++) 
		r_i_a[i] = (long double)(atof(argv[1+3+i]));
	
#else
#ifndef NEW_BASIS
	#define NU 20
	int nu = NU;
	long double r_i_a[NU] = 	 {.2, .4, .65, .8, .9, 
				1.5, 1.7, 
				2.1, 2.4, 2.7,
				3.2, 3.5, 3.8, 
				4.2, 4.5, 4.9, 5.2, 5.8, 6.3, 7.0};
#endif
#endif
#elif SYSTEM == 2	// Dipolare repulsive harte Scheibchen
	out ("\n-> Dipolare harte Scheibchen");

#ifndef BASIS_LINE
	if ((int)argc != 1+5) 
		die ("\n<cmd> <density> <temperature> <b> <rmax> <m>\n");
#else
	if ((int)argc <= 1+5)
		die ("\n<cmd> <density> <temperature> <b> <rmax> <m> [p_...] \n");
#endif
	p = (long double) atof (argv[1]);
	temperature = (long double) atof (argv[2]);
	b_field = (long double) atof (argv[3]) * powl (10.0, -3.0);

	gitter = (unsigned int)(atoi(argv[5]));//1.5*1024;//1.5*1024;
	Rmax = (unsigned int)(atoi(argv[4]));//30.0;//50.0;

	out (" p:%Lf temperature:%Lf b:%Lf Rmax:%Lf gitter:%d", p, temperature, b_field, Rmax, gitter);
#ifdef BASIS_LINE
	int NU = (int)argc - (1+5);
	int nu = NU;
	long double *r_i_a = (long double*) alloc_matrix_1 ((unsigned int)NU);
	for (i = 0; i < NU; i++) 
		r_i_a[i] = (long double)(atof(argv[6+i]));
#else
#ifndef NEW_BASIS
	#define NU 29// 34 
	
	int nu = NU;
	long double r_i_a[NU] ={.2, .4, .65, .9, 1.5, 1.7, 2.1, 2.3, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25.0};
/*	 {.2, .35, .5, .65, .8, .9, 
				1.2, 1.35, 1.5, 1.65, 1.8, 1.95,
				2.2, 2.35, 2.5, 2.65, 2.8, 2.95,
				3.2, 3.35, 3.5, 3.65, 3.8, 3.95, 
				4.2, 4.5, 4.9, 5.2, 5.8, 6.3, 7.0, 8.0, 9.0, 12.5};
*/				
#endif
#endif
#elif SYSTEM == 3	// LJ
	out ("\n-> LJ ");
	if ((int) argc != 1+2)
		die ("\n<cmd> <density> <temperature>\n");
	p = (long double) atof (argv[1]);
	temperature = (long double) atof (argv[2]);
	b_field = 0.0;

	gitter = .25*1024;
	Rmax = 10.0;

	out (" p:%Lf temperature:%Lf Rmax:%Lf gitter:%d", p, temperature, Rmax, gitter);
#ifndef NEW_BASIS
	#define NU 20 //13
	int nu = NU;
	long double r_i_a[NU] = {.05, .1, .15, .2, .3, .4, .7, .9, 1.2, 1.5, 2.0, 2.3, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 7.0, 9.0 };
		// {.2, .4, .65, .95, 1.2, 1.5, 1.7, 2.0, 2.5, 3.0, 4.5, 6.7, 10.0};
#endif
#elif SYSTEM == 4	// Dipolare attraktive harte Scheibchen
	out ("\n-> Dipolare attraktive harte Scheibchen");

	if ((int)argc != 1+3) 
		die ("\n<cmd> <density> <temperature> <b>\n");
	p = (long double) atof (argv[1]);
	temperature = (long double) atof (argv[2]);
	b_field = (long double) atof (argv[3]) * powl (10.0, -3.0);

	gitter = 2*1024;
	Rmax = 30.0;

	out (" p:%Lf temperature:%Lf b:%Lf Rmax:%Lf gitter:%d", p, temperature, b_field, Rmax, gitter);
#ifndef NEW_BASIS
	#define NU 28
	int nu = NU;
	long double r_i_a[NU] = {.1, .2, .3, .5, .7, .9, 1.2, 1.5, 1.7, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0, 4.3, 4.7, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0};
#endif
#endif
	
	
	// Orte der Basisfkt. festlegen
#ifdef NEW_BASIS
	out ("\n-> Neue Methode der Basisfuntionserzeugung");
	int nu = 17;
	long double *r_i_a = malloc ((size_t)(nu*sizeof(long double)));
	long double *r_i_a_reservoir = malloc ((size_t)((nu*(nu+1)/2)*sizeof(long double))),
	     *r_i_a_ptr = r_i_a_reservoir;
	// nu = 17
	r_i_a_reservoir[0] = .1;
	r_i_a_reservoir[1] = .2;
	r_i_a_reservoir[2] = .3;
	r_i_a_reservoir[3] = .45;
	r_i_a_reservoir[4] = .65;
	r_i_a_reservoir[5] = .9;
	r_i_a_reservoir[6] = 1.1;
	r_i_a_reservoir[7] = 1.5;
	r_i_a_reservoir[8] = 1.9;
	r_i_a_reservoir[9] = 2.2;
	r_i_a_reservoir[10] = 2.5;
	r_i_a_reservoir[11] = 2.7;
	r_i_a_reservoir[12] = 3.2;
	r_i_a_reservoir[13] = 4.5;
	r_i_a_reservoir[14] = 6.0;
	r_i_a_reservoir[15] = 8.0;
	r_i_a_reservoir[16] = 10.0;
	// nu = 16
	r_i_a_reservoir[17+0] = .1;
	r_i_a_reservoir[17+1] = .2;
	r_i_a_reservoir[17+2] = .3;
	r_i_a_reservoir[17+3] = .45;
	r_i_a_reservoir[17+4] = .65;
	r_i_a_reservoir[17+5] = .9;
	r_i_a_reservoir[17+6] = 1.1;
	r_i_a_reservoir[17+7] = 1.5;
	r_i_a_reservoir[17+8] = 1.9;
	r_i_a_reservoir[17+9] = 2.5;
	r_i_a_reservoir[17+10] = 3.0;
	r_i_a_reservoir[17+11] = 3.5;
	r_i_a_reservoir[17+12] = 4.5;
	r_i_a_reservoir[17+13] = 6.5;
	r_i_a_reservoir[17+14] = 8.0;
	r_i_a_reservoir[17+15] = 10.0;
	// nu = 15
	r_i_a_reservoir[17+16+0] = .1;
	r_i_a_reservoir[17+16+1] = .2;
	r_i_a_reservoir[17+16+2] = .3;
	r_i_a_reservoir[17+16+3] = .45;
	r_i_a_reservoir[17+16+4] = .65;
	r_i_a_reservoir[17+16+5] = .9;
	r_i_a_reservoir[17+16+6] = 1.1;
	r_i_a_reservoir[17+16+7] = 1.5;
	r_i_a_reservoir[17+16+8] = 1.9;
	r_i_a_reservoir[17+16+9] = 2.5;
	r_i_a_reservoir[17+16+10] = 3.0;
	r_i_a_reservoir[17+16+11] = 3.5;
	r_i_a_reservoir[17+16+12] = 4.5;
	r_i_a_reservoir[17+16+13] = 6.5;
	r_i_a_reservoir[17+16+14] = 8.5;
	// nu = 14
	r_i_a_reservoir[17+16+15+0] = .1;
	r_i_a_reservoir[17+16+15+1] = .2;
	r_i_a_reservoir[17+16+15+2] = .3;
	r_i_a_reservoir[17+16+15+3] = .45;
	r_i_a_reservoir[17+16+15+4] = .65;
	r_i_a_reservoir[17+16+15+5] = .9;
	r_i_a_reservoir[17+16+15+6] = 1.1;
	r_i_a_reservoir[17+16+15+7] = 1.5;
	r_i_a_reservoir[17+16+15+8] = 1.9;
	r_i_a_reservoir[17+16+15+9] = 2.5;
	r_i_a_reservoir[17+16+15+10] = 3.0;
	r_i_a_reservoir[17+16+15+11] = 4.5;
	r_i_a_reservoir[17+16+15+12] = 5.5;
	r_i_a_reservoir[17+16+15+13] = 8.5;
	// nu = 13
	r_i_a_reservoir[17+16+15+14+0] = .1;
	r_i_a_reservoir[17+16+15+14+1] = .2;
	r_i_a_reservoir[17+16+15+14+2] = .3;
	r_i_a_reservoir[17+16+15+14+3] = .45;
	r_i_a_reservoir[17+16+15+14+4] = .65;
	r_i_a_reservoir[17+16+15+14+5] = .9;
	r_i_a_reservoir[17+16+15+14+6] = 1.1;
	r_i_a_reservoir[17+16+15+14+7] = 1.5;
	r_i_a_reservoir[17+16+15+14+8] = 1.9;
	r_i_a_reservoir[17+16+15+14+9] = 2.5;
	r_i_a_reservoir[17+16+15+14+10] = 3.0;
	r_i_a_reservoir[17+16+15+14+11] = 5.5;
	r_i_a_reservoir[17+16+15+14+12] = 8.5;
	// nu = 12
	r_i_a_reservoir[17+16+15+14+13+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+1] = .2;
	r_i_a_reservoir[17+16+15+14+13+2] = .3;
	r_i_a_reservoir[17+16+15+14+13+3] = .45;
	r_i_a_reservoir[17+16+15+14+13+4] = .65;
	r_i_a_reservoir[17+16+15+14+13+5] = .9;
	r_i_a_reservoir[17+16+15+14+13+6] = 1.1;
	r_i_a_reservoir[17+16+15+14+13+7] = 1.5;
	r_i_a_reservoir[17+16+15+14+13+8] = 1.9;
	r_i_a_reservoir[17+16+15+14+13+9] = 2.5;
	r_i_a_reservoir[17+16+15+14+13+10] = 4.7;
	r_i_a_reservoir[17+16+15+14+13+11] = 8.5;
	// nu = 11
	r_i_a_reservoir[17+16+15+14+13+12+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+12+1] = .25;
	r_i_a_reservoir[17+16+15+14+13+12+2] = .3;
	r_i_a_reservoir[17+16+15+14+13+12+3] = .45;
	r_i_a_reservoir[17+16+15+14+13+12+4] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+5] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+6] = 1.1;
	r_i_a_reservoir[17+16+15+14+13+12+7] = 1.5;
	r_i_a_reservoir[17+16+15+14+13+12+8] = 2.5;
	r_i_a_reservoir[17+16+15+14+13+12+9] = 4.5;
	r_i_a_reservoir[17+16+15+14+13+12+10] = 8.0;
	// nu = 10
	r_i_a_reservoir[17+16+15+14+13+12+11+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+12+11+1] = .25;
	r_i_a_reservoir[17+16+15+14+13+12+11+2] = .45;
	r_i_a_reservoir[17+16+15+14+13+12+11+3] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+11+4] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+5] = 1.4;
	r_i_a_reservoir[17+16+15+14+13+12+11+6] = 2.2;
	r_i_a_reservoir[17+16+15+14+13+12+11+7] = 3.0;
	r_i_a_reservoir[17+16+15+14+13+12+11+8] = 4.5;
	r_i_a_reservoir[17+16+15+14+13+12+11+9] = 8.0;
	// nu = 9
	r_i_a_reservoir[17+16+15+14+13+12+11+10+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+1] = .25;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+2] = .45;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+3] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+4] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+5] = 1.2;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+6] = 1.9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+7] = 2.5;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+8] = 4.9;
	// nu = 8
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+1] = .3;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+2] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+3] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+4] = 1.2;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+5] = 1.7;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+6] = 2.5;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+7] = 4.9;
	// nu = 7
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+0] = .1;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+1] = .3;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+2] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+3] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+4] = 1.2;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+5] = 1.9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+6] = 3.5;
	// nu = 6
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+0] = .25;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+1] = .65;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+2] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+3] = 1.2;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+4] = 2.4;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+5] = 3.0;
	// nu = 5
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+6+0] = .5;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+6+1] = .9;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+6+2] = 1.5;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+6+3] = 2.3;
	r_i_a_reservoir[17+16+15+14+13+12+11+10+9+8+7+6+4] = 3.0;

	Rmax = 100.0;
#endif

	out ("\n-> Basen (%d):", nu);
	for (i = 0; i < nu; i++)
		out (" %.1Lf", r_i_a[i]);

	out ("\n-> Mittlerer Abstand: %lf", pow(p,-.5));
	if (pow(p,-.5) < 1.0)
		die ("!\n");


	/**********************************************
	 * Funktionen für Iteration und Zerlegung initialisieren
	 **********************************************/
	out ("\n-> Funktionen initialisieren");
	long double *C_r = alloc_matrix_1 (gitter);
	long double *C_q = alloc_matrix_1 (gitter);
	long double *Gamma_r = alloc_matrix_1 (gitter);
	long double *Gamma_q = alloc_matrix_1 (gitter);

	
	/******************
	 * Besselfunktionen
	 ******************/
	out ("\n-> Besselfunktionen");
	// Punkteraster: nicht äquidistant (siehe Lado)
	long double Qmax = gsl_sf_bessel_zero_J0 (gitter+1) / Rmax;
	long double *r_ = alloc_matrix_1 (gitter);
	for (i = 0; i < gitter; i++) 
		r_[i] = (long double) (gsl_sf_bessel_zero_J0 (i+1) / Qmax);
	long double *q_ = alloc_matrix_1 (gitter);
	for (i = 0; i < gitter; i++) 
		q_[i] = (long double) (gsl_sf_bessel_zero_J0 (i+1) / Rmax);
	
	// Struktur für FBT 1
	t_fbt fbt1 = {
		.gitter = gitter,
		.cons = 4.0*M_PIl/(Qmax*Qmax),
		.bessel = alloc_matrix_2 (gitter, gitter)
	};
	for (m = 0; m < gitter; m++) for (i = 0; i < gitter; i++) 
		fbt1.bessel[m][i] = j0l (q_[m] * r_[i]) / powl (j1l (r_[i] * Qmax), 2.0);

	// Struktur für FBT 2
	t_fbt fbt2 = {
		.gitter = gitter,
		.cons = 1.0/(M_PIl*Rmax*Rmax),
		.bessel = alloc_matrix_2 (gitter, gitter)
	};
	for (m = 0; m < gitter; m++) for (i = 0; i < gitter; i++) 
		fbt2.bessel[m][i] = j0l (q_[m] * r_[i]) / powl (j1l (q_[m] * Rmax), 2.0);

	
	/***************************
	 * Mayerfunktionen berechnen
	 ***************************/
	out ("\n-> Mayerfunktion");
	long double *mayer = alloc_matrix_1 (gitter);
	long double *Vr = alloc_matrix_1 (gitter); 
	return_mayer (mayer, r_, Vr, gitter, p, temperature, b_field);
	write_function (mayer, gitter, r_, "mayer", 0); 
#ifdef _MAYER_
	goto clean;
#endif


	/************
	 * Startwerte
	 ************/
	out ("\n-> Startwerte");
#if STARTVALUE == 1
	out (" Gamma(r) == 0");
	for (i = 0; i < gitter; i++) 
		Gamma_r[i] = 0.0; 
#elif STARTVALUE == 2
	out (" Gamma(r) = 'Fitfunktion'");
	// In reduzierten Einheiten für harte Scheibchen einer Sorte bei red. p = 0.7 gilt gehnähert:
	// Gamma(r) = 20*exp(-x^2/4)
	/* Fixme:
	 * Diese Funktion stimmt in den red. Einheiten nicht mehr
	 */
	for (i = 0; i < gitter; i++) 
		Gamma_r[i] = 20.0*expl (-powl (r_[i], 2.0)/4.0); 
#elif STARTVALUE == 3
	out (" Gamma(r) einlesen aus 'startvalue'");
	read_function (Gamma_r, r_, gitter, "startvalue");
#elif STARTVALUE == 4
	out (" Gamma(r) == Mayerfkt");
	for (i = 0; i < gitter; i++)
		Gamma_r[i] = mayer[i];
#endif


#ifdef DYN_BASIS
	r_i_a_ptr = r_i_a_reservoir;
	for (n = 0; n < nu; n++, r_i_a_ptr++)
		r_i_a[n] = *r_i_a_ptr;
	out ("\n-> %d Basen: ", nu);
	for (n = 0; n < nu; n++)
		out (" %.3Lf", r_i_a[n]);
#endif

		
	/*****************
	 * Zerlegungsbasis 
	 *****************/
#ifdef NEW_BASIS
	long double **P_ = alloc_matrix_2 (nu, gitter),	// Cleanup: Das könnte man in externe routine machen 
		**Q_ = alloc_matrix_2 (nu, gitter);
	pq_init (r_i_a, (unsigned int)nu, r_, (unsigned int)gitter, P_, Q_);
#else
	/* Cleanup: In Fkt- pq_gen */
	out ("\n-> Zerlegungsbasis ");
	// Stützpunkte wählen
#ifndef BASIS_LINE
	int i_[NU];   // Die nu i_alphas aus den Maximalstellen bestimmen
#else
	int *i_ = (int *) calloc ((size_t)NU, (size_t)sizeof (int));
#endif
	for (n = 0; n < nu; n++) 
		i_[n] = find_i (r_, r_i_a[n], gitter);
	
	// *** P_n^i konstruieren ***
	out ("P_n^i ");
	long double **P_ = alloc_matrix_2 (nu, gitter);
	for (n = 0; n < nu; n++) 
		if (n == 0) {   // 1. Basis
			for (i = 0; i < gitter; i++) {
				P_[0][i] = 0.0;
				//if (i <= i_[0]) P_[0][i] = (long double) (i_[0]-i) / (i_[0]);
				if (i <= i_[1]) P_[0][i] = (long double) (i_[1]-i) / (long double)(i_[1]);
			}
		} else if (n == 1) {	// 2. Basis
			for (i = 0; i < gitter; i++) {
				P_[n][i] = 0.0;
				if (0 <= i && i <= i_[1]) P_[n][i] = (long double)(i)/(long double)(i_[1]);
				else if (i_[1] <= i && i <= i_[2]) P_[n][i] = (long double)(i_[2]-i)/(long double)(i_[2]-i_[1]);
			}
	 	} else {   // Alle anderen Basen
			for (i = 0; i < gitter; i++) {
				P_[n][i] = 0.0;
				if (i_[n-1] <= i && i <= i_[n]) P_[n][i] = fabsl ((long double) (i-i_[n-1]) / (i_[n-1]-i_[n]));
				if (i_[n] <= i && i <= i_[n+1])	P_[n][i] = fabsl ((long double) (i_[n+1]-i) / (i_[n+1]-i_[n]));
			}
		}

	// *** Q_n^i konstruieren ***
	out ("Q_n^i ");
	long double **Q_ = alloc_matrix_2 (nu, gitter);
	
	// Dazu die Matrix (PP)_mn = sum_i(P_m^i*P_n^i) berechnen
	long double sum;
	long double **PP = alloc_matrix_2 (nu, nu);
	for (m = 0; m < nu; m++) for (n = 0; n < nu; n++) {
		sum = 0.0;
		for (i = 0; i < gitter; i++)
			sum += P_[m][i]*P_[n][i];
		PP[m][n] = sum;
	}

	// Dazu die Inserve RR von PP berechnen
	long double **RR = alloc_matrix_2 (nu, nu);
	invert (PP, RR, nu);

	// Und schließlich die Qs berechnen
	for (m = 0; m < nu; m++) 
		for (i = 0; i < gitter; i ++) {
			sum = 0.0;
			for (n = 0; n < nu; n++) 
				sum += RR[m][n]*P_[n][i];
			Q_[m][i] = sum;
		}
#endif

#define _BASIS_
#ifdef _BASIS_
	FILE *basis_file = fopen ("basis", "w");
	if (!basis_file)
		die ("\nbasis_file: basis\n");
	for (i = 0; i < gitter; i++) {
		fprintf (basis_file, "\n%Lf ", r_[i]);
		for (m = 0; m < nu; m++) 
			fprintf (basis_file, "%Lf ", P_[m][i]);
		for (m = 0; m < nu; m++)
			fprintf (basis_file, "%Lf ", Q_[m][i]);
	}
	fclose (basis_file);
#endif

	// Strukturen für Zerlegung / Komposition erstellen
	t_composite comp = {
		.function = Gamma_r,
		.coarse = alloc_matrix_1 (nu),
		.fine = alloc_matrix_1 (gitter),
		.P_ = P_,
		.Q_ = Q_,
		.gitter = gitter,
		.nu = nu
	};

	t_composite decomp = {
		.function = Gamma_r,
		.coarse = alloc_matrix_1 (nu),
		.fine = alloc_matrix_1 (gitter),
		.P_ = P_,
		.Q_ = Q_,
		.gitter = gitter,
		.nu = nu
	};
	
	// Aufräumen
#ifndef NEW_BASIS
	free_matrix_2 (RR, nu);
	free_matrix_2 (PP, nu);
#endif


	
	/***************
	 * Hauptschleife
	 ***************/
	out ("\n-> Hauptschleife");
	int pic = 0;
	long double norm_fine = 10.0, norm_fine_old = 10.0;
	long double norm_fine_max = 1e-10;	// Konvergenzkriterium
	long double norm_sum = 0.0;

	// Ausgabe vorbereiten
	t_functions functions = {
		.Gamma_r = Gamma_r,
		.Gamma_q = Gamma_q,
		.C_r = C_r,
		.C_q = C_q,
		.r_ = r_,
		.q_ = q_,
		.gitter = &gitter,
		.iteration = &pic,
		.p = &p
	};

	// Startwerte ausgeben
#define __STARTFUNCTION
#ifdef __STARTFUNCTION
	function_out (functions);
#endif

	// Ausgabe von Grob- und Feinanteil vorbereiten
#ifdef _cf_
	t_coarse_fine coarse_fine = {
		.coarse = decomp.coarse,
		.fine = decomp.fine,
		.coarse_dim = &nu,
		.fine_dim = gitter,
		.nr = &nr,
		.pic = &pic,
		.r_ = r_,
		.P_ = P_
	};
	coarse_fine.coarse_dim = coarse_fine.coarse_dim;	// Gegen Warnung
#endif

	
	/**************
	 * Testroutinen
	 **************/
#ifdef TEST
	// Ende der Testroutinen
	goto clean;
#endif

	// Erste Zerlegung
	out ("\n-> Erste Zerlegung");
	decomposite (comp);
#ifdef DYN_BASIS_1
	for (i = 0; i < gitter; i++)
		decomp.fine[i] = comp.fine[i];
	for (n = 0; n < nu; n++)
		decomp.coarse[n] = comp.coarse[n];
#endif

	// LM-Date
	lm_data_type data = {	/* Fixme: Nur diese eine Struktur verwenden */
		.Cr = C_r, .Cq = C_q, .Gammar = Gamma_r, .Gammaq = Gamma_q,
		.mayer = mayer, .Vr = Vr,
		.r_ = r_,
		.p = p,
		.gitter = gitter, .nu = nu,
		.fbt1 = fbt1, .fbt2 = fbt2,
		.comp = {
			.function = Gamma_r,
			.coarse = alloc_matrix_1 (nu),
			.fine = alloc_matrix_1 (gitter),
			.P_ = P_, .Q_ = Q_,
			.gitter = gitter, .nu = nu
		},
		.decomp = {			
			.function = Gamma_r,
			.coarse = alloc_matrix_1 (nu),
			.fine = alloc_matrix_1 (gitter),
			.P_ = P_, .Q_ = Q_,
			.gitter = gitter, .nu = nu
		}
	};

	// Beginn:	
	time_t start, now;
	time (&start);

	j = v = 1.0; 	// Gegen Warnung
	
	lm_control_type control;
	lm_initialize_control (&control);

	//function_out (functions);
	// Startwerte für LM Genauigkeit:
	control.epsilon = control.ftol = control.xtol = control.gtol = 1e-25;

	
	// Grobanteilgenauigkeit dynamisch anpassen
#ifdef DYN_COARSE
	out ("\n-> Grobanteilgenauigkeit adaptiv anpassen");
	// Schrittweitenparameter
	double dyn_coarse_step = 1e-1,			// Multiplikator für adaptive Genauigkeit im LM
		dyn_coarse = 1e-3,			// Startgenauigkeit für LM
		dyn_coarse_max = 1e-12;		// Maximale Genauigkeit für LM
	double norm_fine_max_dyn_coarse = 1e-10;		// Startkriterium für Feinanteilgenauigkeit
	// Startwerte
	control.epsilon = dyn_coarse;
	control.ftol = dyn_coarse;
	control.xtol = dyn_coarse;
	control.gtol = dyn_coarse;
#endif
#ifdef DYN_COARSE_2
	out ("\n-> Grobanteilgenauigkeit adaptiv anpassen (2)");
	// Schrittweitenparameter
	double dyn_coarse_step = 1e-1,
		dyn_coarse = 1e-12,
		dyn_coarse_min = 1e-1;
	double norm_fine_max_dyn_coarse = 1e-1;
	// Startwerte
	control.epsilon = dyn_coarse;
	control.ftol = dyn_coarse;
	control.xtol = dyn_coarse;
	control.gtol = dyn_coarse;
#endif
#ifdef DYN_BASIS
	double norm_fine_max_dyn_basis = 1e-5;
#endif
#ifdef DYN_BASIS_1
	double norm_fine_max_dyn_basis = 1e-3;
#endif

#ifdef DYN_BASIS_1
	out ("\n-> Neu Basis (1)");
			// Gamma zusammensetzen
			composite (data.decomp);

			// Alte Basen & Fkt. löschen
			free_matrix_2 (P_, nu);
			free_matrix_2 (Q_, nu);

			// Neues nu & Grobanteile
			nu=7;
			data.nu = data.comp.nu = data.decomp.nu = comp.nu = decomp.nu = nu; 	// Fixme: Das hier geht mit Zeigern...
			data.comp.coarse = (long double *) realloc (data.comp.coarse, (size_t)(nu*sizeof(long double)));
			data.decomp.coarse = (long double *) realloc (data.decomp.coarse, (size_t)(nu*sizeof(long double)));
			comp.coarse = (long double *) realloc (comp.coarse, (size_t)(nu*sizeof(long double)));
			decomp.coarse = (long double *) realloc (decomp.coarse, (size_t)(nu*sizeof(long double)));
			r_i_a = (long double *)malloc ((size_t)(nu*sizeof(long double)));	/* Fixme: geht das einfacher? */
			// Neue Basis
			for (m = 17, j = 0; m > nu; m--)
				j += m;
			for (n = 0, r_i_a_ptr = r_i_a_reservoir+j; n < nu; n++, r_i_a_ptr++)
				r_i_a[n] = *r_i_a_ptr;

			// Neue Basen
			P_ = alloc_matrix_2 (nu, gitter);
			Q_ = alloc_matrix_2 (nu, gitter);
			pq_init (r_i_a, (unsigned int)nu, r_, (unsigned int)gitter, P_, Q_);
			data.comp.P_ = data.decomp.P_ = comp.P_ = decomp.P_ = P_;
			data.comp.Q_ = data.decomp.Q_ = comp.Q_ = decomp.Q_ = Q_;

			// Zerlegung in neue Basen
			decomposite (data.decomp);
			for (i = 0; i < gitter; i++)
				comp.fine[i] = data.decomp.fine[i];

			// Ausgabe
			out ("\n-> %d Basen:", nu);
			for (n = 0; n < nu; n++)
				out (" %.3Lf", r_i_a[n]);
#endif



	// Hauptschleife
	out ("\n-> Hauptschleife");
	unsigned int pic_max = 10000;
	unsigned int out_period = 1;	// Periode zum Rausschreiben
	double out_period_norm = 1e-1;
	do {
		pic++;
		out ("\n%d", pic);

		// Grob
		for (i = 0; i < gitter; i++)
			data.comp.fine[i] = comp.fine[i];	// Feinanteil für LM; Grobanteil wird ja eben variiert
	
		out (" lm %d ", nu);
		lm_minimize (nu, nu, (double *)comp.coarse, lm_evaluate_default, lm_print_default, &data, &control);

		norm_fine_old = norm_fine;
		norm_fine = norm_sum = 0.0;
		for (i = 0; i < gitter; i++) {
			norm_fine += powl ((comp.fine[i]-data.decomp.fine[i]), 2.0);
			norm_sum += data.decomp.fine[i];
	//		comp.fine[i] = data.decomp.fine[i];
		}
		norm_fine = sqrt (norm_fine)/norm_sum;
		out (" fine: %g", (double)norm_fine);
		
		if (norm_fine_old < norm_fine) {
			out (" !! ");
			for (i = 0; i < gitter; i++) 
				comp.fine[i] = .1*data.decomp.fine[i] + (1.0-.1)*comp.fine[i];
			//norm_fine = norm_fine_old;
		} else {
			for (i = 0; i < gitter; i++)
				comp.fine[i] = data.decomp.fine[i];
		}

#ifdef DYN_BASIS_1
		if (norm_fine < norm_fine_max_dyn_basis && nu < 16) {
			// Gamma zusammensetzen
			composite (data.decomp);

			// Alte Basen & Fkt. löschen
			free_matrix_2 (P_, nu);
			free_matrix_2 (Q_, nu);

			// Neues nu & Grobanteile
			nu++;
			data.nu = data.comp.nu = data.decomp.nu = comp.nu = decomp.nu = nu; 	// Fixme: Das hier geht mit Zeigern...
			data.comp.coarse = (long double *) realloc (data.comp.coarse, (size_t)(nu*sizeof(long double)));
			data.decomp.coarse = (long double *) realloc (data.decomp.coarse, (size_t)(nu*sizeof(long double)));
			comp.coarse = (long double *) realloc (comp.coarse, (size_t)(nu*sizeof(long double)));
			decomp.coarse = (long double *) realloc (decomp.coarse, (size_t)(nu*sizeof(long double)));
			r_i_a = (long double *)malloc ((size_t)(nu*sizeof(long double)));	/* Fixme: geht das einfacher? */
			// Neue Basis
			for (m = 17, i = 0; m > nu; m--)
				i += m;
			for (n = 0, r_i_a_ptr = r_i_a_reservoir+i; n < nu; n++, r_i_a_ptr++)
				r_i_a[n] = *r_i_a_ptr;

			// Neue Basen
			P_ = alloc_matrix_2 (nu, gitter);
			Q_ = alloc_matrix_2 (nu, gitter);
			pq_init (r_i_a, (unsigned int)nu, r_, (unsigned int)gitter, P_, Q_);
			data.comp.P_ = data.decomp.P_ = comp.P_ = decomp.P_ = P_;
			data.comp.Q_ = data.decomp.Q_ = comp.Q_ = decomp.Q_ = Q_;

			// Zerlegung in neue Basen
			decomposite (data.decomp);
			for (i = 0; i < gitter; i++)
				comp.fine[i] = data.decomp.fine[i];

			// Ausgabe
			out ("\n-> %d Basen:", nu);
			for (n = 0; n < nu; n++)
				out (" %.3Lf", r_i_a[n]);
		}
#endif

#ifdef DYN_BASIS
		if (norm_fine < norm_fine_max_dyn_basis && nu == 14) {
			// Gamma zusammensetzen
			composite (data.decomp);

			// Alte Basen & Fkt. löschen
			free_matrix_2 (P_, nu);
			free_matrix_2 (Q_, nu);

			// Neues nu & Grobanteile
			nu = 14;
			data.nu = data.comp.nu = data.decomp.nu = comp.nu = decomp.nu = nu; 	// Fixme: Das hier geht mit Zeigern...
			data.comp.coarse = (long double *) realloc (data.comp.coarse, (size_t)(nu*sizeof(long double)));
			data.decomp.coarse = (long double *) realloc (data.decomp.coarse, (size_t)(nu*sizeof(long double)));
			comp.coarse = (long double *) realloc (comp.coarse, (size_t)(nu*sizeof(long double)));
			decomp.coarse = (long double *) realloc (decomp.coarse, (size_t)(nu*sizeof(long double)));
			r_i_a = (long double *)malloc ((size_t)(nu*sizeof(long double)));	/* Fixme: geht das einfacher? */
			for (n = 0; n < nu; n++, r_i_a_ptr++)
				r_i_a[n] = *r_i_a_ptr;

			// Neue Basen
			P_ = alloc_matrix_2 (nu, gitter);
			Q_ = alloc_matrix_2 (nu, gitter);
			pq_init (r_i_a, (unsigned int)nu, r_, (unsigned int)gitter, P_, Q_);
			data.comp.P_ = data.decomp.P_ = comp.P_ = decomp.P_ = P_;
			data.comp.Q_ = data.decomp.Q_ = comp.Q_ = decomp.Q_ = Q_;

			// Zerlegung in neue Basen
			decomposite (data.decomp);
			for (i = 0; i < gitter; i++)
				comp.fine[i] = data.decomp.fine[i];

			// Ausgabe
			out ("\n-> %d Basen:", nu);
			for (n = 0; n < nu; n++)
				out (" %.3Lf", r_i_a[n]);
		}
#endif
#ifdef DYN_COARSE
		if (norm_fine < dyn_coarse && dyn_coarse > dyn_coarse_max) { //norm_fine_max_dyn_coarse && dyn_coarse > dyn_coarse_max) {
			dyn_coarse *= dyn_coarse_step;
			control.epsilon = dyn_coarse;
			control.ftol = dyn_coarse;
			control.xtol = dyn_coarse;
			control.gtol = dyn_coarse;
			out ("\n-> Grobanteilgenauigkeit anpassen %g\n", dyn_coarse);
#ifdef DEBUG
//			goto clean;
#endif
		}
#endif
#ifdef DYN_COARSE_2
		if (norm_fine < norm_fine_max_dyn_coarse && dyn_coarse < dyn_coarse_min) {
			dyn_coarse /= dyn_coarse_step;
			control.epsilon = dyn_coarse;
			control.ftol = dyn_coarse;
			control.xtol = dyn_coarse;
			control.gtol = dyn_coarse;
			out ("\n-> Grobanteilgenauigkeit anpassen %.15lf", dyn_coarse);
		}
#endif

#ifdef DEBUG
		function_out (functions);
#else
		if ((pic % 100) == 0) 
			function_out (functions);
#endif


		/* Ausgabe */
		if (norm_fine < out_period_norm) {
			function_out (functions);
			out_period_norm = .1*norm_fine;
			out_period *= 10;
		} 
		if ((pic % out_period) == 0) 
			function_out (functions);

		time (&now);
		out (" time: %.0lf", difftime (now, start));
#ifdef DYN_COARSE
	} while (norm_fine > norm_fine_max && pic < pic_max);// && dyn_coarse > dyn_coarse_max);
#else
	} while (norm_fine > norm_fine_max && pic < pic_max);
#endif

	// Cleanup
#ifdef DYN_COARSE
	norm_fine_max_dyn_coarse = 0.0; // Gegen Warnung
#endif
	goto clean;	// Gegen Warnung
clean:
	out ("\n-> Cleanup");
	function_out (functions);

	/***********
	 * Aufräumen
	 ***********/
	free_matrix_1 (C_r);
	free_matrix_1 (C_q);
	free_matrix_1 (mayer);
	free_matrix_1 (Vr);
	free_matrix_1 (Gamma_r);
	free_matrix_1 (Gamma_q);
#ifdef PIC_AFTER
	free_matrix_1 (Gamma__);
#endif

	free_matrix_1 (r_);
	free_matrix_1 (q_);
	free_matrix_2 (fbt1.bessel, gitter);
	free_matrix_2 (fbt2.bessel, gitter);

#ifdef NEW_BASIS
	free_matrix_1 (r_i_a);
	free_matrix_1 (r_i_a_reservoir);
#endif

	free_matrix_2 (P_, nu);
	free_matrix_2 (Q_, nu);

	free_matrix_1 (data.comp.fine);
	free_matrix_1 (data.comp.coarse);
	free_matrix_1 (data.decomp.fine);
	free_matrix_1 (data.decomp.coarse);
	
	free_matrix_1 (comp.fine);
	free_matrix_1 (comp.coarse);
	free_matrix_1 (decomp.fine);
	free_matrix_1 (decomp.coarse);

#ifdef BASIS_LINE
	free (i_);
	free (r_i_a);
#endif
	
	out ("\n-> Fertig\n");
	if (pic >= pic_max)
		return (-1);

	return (pic);
}
