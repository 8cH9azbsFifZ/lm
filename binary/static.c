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
#include <stdarg.h>	// Für out / die
#define __USE_GNU 1 // für M_PIl: long double PI aus libmath
#define _GNU_SOURCE 1 // Für j0l: long double Bessel aus libmath
double j0(double x);	// Gibt sonst Warnungen beim Compilieren
double j1(double y);
#include <math.h>
#include <getopt.h>
#include <string.h>

#include <time.h>

#include <assert.h>
#define NDEBUG	// Für assert

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


/******************************
 * Schalter für den Algorithmus
 ******************************/
/* Fixme: Potentiale trennen */
/* System:
 * 1: Harte Scheibchen
 * 2: Mix harter Scheibchen 
 * 3: Dipolare WW
 * 4: LJ
 * 5: exponentieller Abfall
 * 6: auch die suszeptibilitäten variieren
 * 7: Dipolar und attraktiv
 */
#define SYSTEM 3
//#define SYSTEM_3	// Debug
//#define FINE_RESCALE

/* Approximation
 * 1: PY
 * 2: HNC
 */
#define APPROXIMATION 1

#define BASIS_LINE
#define PACKING_FRACTION 	// eta angeben

//#define DYN_COARSE   // Grobanteilwachstum beschränken
//#define SYM 2    // Gamma(r) nach OZ immer wieder symmetrisieren 
		//1: Alle 4 Elemente gleich machen
		//2: Symmetrisieren
//#define ASYMPTOTE   // Gamma(r) = -1-C(r)   r<1
#define STARTVALUE 2   // Startwert:   
			// 1: Gamma(r)=Startfkt 
			// 2: Gamma(r)=0 
			// 3: h(r) einlesen 
			// 4: Gamma(r)=.1
			// 5: Gamma(r) einlesen
#define SYMMETRIC_FT

/* Testroutinen:
 * 1 = FBT
 * 2 = Zerlegung
 * 3 = Startgamma mit cq_to_sq
 * 5 = Invertierung
 * 6 = Zerlegung wie im echten Leben...
 * 7 = Startgamma einlesen
 */
//#define DEBUG 100   // # Iterationen
//#define TEST 7

#define MARQUARDT

//#define DYN_COARSE

//#define T_PATH




/*******
 * Typen
 *******/
/* Struktur für eine Matrix in einem 1d Array mit Zeigern für diese */
typedef struct {
	double *matrix;		// Start des Arrays
	double *ptr, *ptr1;	// Zeiger in dem Array
	unsigned int size;	// Größe des Arrays
} t_matrix;

/* Konstenaten & Fkt. für die FT */
typedef struct {
	t_matrix mr;
	t_matrix mq;
	t_matrix bessel;
	int gitter;	// = sqrt (mr.size)
	double cons;
} t_fbt;

/* Konstanten, Basen & Fkt. für Zerlegung / Komposition */
typedef struct {
	t_matrix function;
	t_matrix coarse, fine;
	t_matrix P_, Q_;
	unsigned int gitter, nu;
	double sum;   // temporäre Zwischensumme
} t_composite;

/* Debugausgabe der Basen */
typedef struct {
	t_matrix P_;
	unsigned int gitter, nu;
	t_matrix r_;
} t_basis_out;

/* Debugausgabe von Grob-/Fein-Anteil */
typedef struct {
	t_matrix coarse;
	t_matrix fine;
	unsigned int *pic;
	unsigned int gitter, nu;
	t_matrix r_;
	t_matrix P_;
} t_coarse_fine;

/* Alle Matrizen und Konstanten für Ausgabe */
typedef struct {
	t_matrix Gammar, Gammaq, Cr, Cq;
	t_matrix r_, q_;
	unsigned int *pic;
	double *x_s, *x_b;
	double *p;
} t_matrices;

/* Wird beim Einlesen einer Matrix zurückgegeben */
typedef struct {
	int rows;
	double stepsize;
} t_read_matrix;

/* Alle Daten für einen Picarsschritt im LM */
typedef struct {
	t_matrix Cr, Cq, Gammar, Gammaq;
	double p_s, p_b;
	t_matrix mayer, Vr;
	t_fbt fbt1, fbt2;
	t_composite comp, decomp;	/* Fixme: In Hauptschleife nur auf Data-Type arbeiten */
} lm_data_type;



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
                        
/* Programm mit Fehlerausgabe beenden */
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
/* Speicherplatz für Matrix allokieren, Zeiger setzen */
t_matrix init_matrix (unsigned int size) {
	t_matrix matrix = {
		.size = size,
		.matrix = (double *) calloc ((size_t)(size), (size_t)(sizeof (double))),
		.ptr = NULL,
		.ptr1 = NULL
	};
	if (!matrix.matrix)
		die ("\nMatrixinit");
	return matrix;
}

/* Speicherplatz einer Matrix freigeben */
void cleanup_matrix (t_matrix matrix) {
	free (matrix.matrix);
}



/*******************
 * Marquardt-Methode
 *******************/
#include "marquardt.c"



/***************************************
 * Routinen für die Picard-Hauptschleife
 ***************************************/
/* PY-Approximation */
void PY (t_matrix Cr, t_matrix Gammar, t_matrix mayer) {
	unsigned int i;
	for (i = 0, Cr.ptr = Cr.matrix, Gammar.ptr = Gammar.matrix, mayer.ptr = mayer.matrix; i < Cr.size; i++, Cr.ptr++, Gammar.ptr++, mayer.ptr++)
		*Cr.ptr = *mayer.ptr*(1.0+*Gammar.ptr);
}

/* HNC-Approximation */
void HNC (t_matrix Cr, t_matrix Gammar, t_matrix Vr) {
	unsigned int i;
	for (i = 0, Cr.ptr = Cr.matrix, Gammar.ptr = Gammar.matrix, Vr.ptr = Vr.matrix; i < Cr.size; i++, Cr.ptr++, Gammar.ptr++, Vr.ptr++)
		*Cr.ptr = exp (*Vr.ptr+*Gammar.ptr)-*Gammar.ptr-1.0;
}

/* FT: Orts- -> Impulsraum */
void FBT1 (t_fbt f) {
//#define _FBT1_
	unsigned int r, q;

	double *ptr11_r, *ptr12_r, *ptr21_r, *ptr22_r;
	double *ptr11_q, *ptr12_q, *ptr21_q, *ptr22_q;

	for (q = 0, f.bessel.ptr = f.bessel.matrix, ptr11_q = f.mq.matrix, ptr12_q = ptr11_q+1, ptr21_q = ptr12_q+1, ptr22_q = ptr21_q+1; 
			q < f.gitter; q++, ptr11_q+=4, ptr12_q+=4, ptr21_q+=4, ptr22_q+=4) {
		*ptr11_q = *ptr12_q = *ptr21_q = *ptr22_q = 0.0;
		for (r = 0, ptr11_r = f.mr.matrix, ptr12_r = ptr11_r+1, ptr21_r = ptr12_r+1, ptr22_r = ptr21_r+1; 
				r < f.gitter; r++, ptr11_r+=4, ptr12_r+=4, ptr21_r+=4, ptr22_r+=4, f.bessel.ptr++) {
			*ptr11_q += *ptr11_r**f.bessel.ptr;
			*ptr12_q += *ptr12_r**f.bessel.ptr;
#ifndef SYMMETRIC_FT
			*ptr21_q += *ptr21_r**f.bessel.ptr;
#endif
			*ptr22_q += *ptr22_r**f.bessel.ptr;
		}
		*ptr11_q *= f.cons;
		*ptr12_q *= f.cons;
#ifdef SYMMETRIC_FT
		*ptr21_q = *ptr12_q;
#else
		*ptr21_q *= f.cons;
#endif
		*ptr22_q *= f.cons;
	}
	/* Fixme: Symmetrie ausnutzen */
}

/* OZ-Gleichung */
void OZ (t_matrix Cq, t_matrix Gammaq, double ps, double pb) {
	unsigned int i;
	double z1;
	Gammaq.ptr = Gammaq.matrix;
	double *ptr11, *ptr12, *ptr21, *ptr22;
	ptr11 = Cq.matrix;
	ptr12 = Cq.matrix+1;
	ptr21 = Cq.matrix+2;
	ptr22 = Cq.matrix+3;
	for (i = 0; i < Gammaq.size; i+=4, Gammaq.ptr+=4, ptr11+=4, ptr12+=4, ptr21+=4, ptr22+=4) {
		z1 = 1.0/(1.0-(*ptr11+*ptr12**ptr21*pb)*ps+*ptr22*pb*(-1.0+*ptr11*ps));
		*(Gammaq.ptr) = z1*(*ptr12**ptr21*pb+*ptr11*(*ptr11+*ptr12**ptr21*pb-*ptr11**ptr22*pb)*ps);
		*(Gammaq.ptr+1) = *ptr12*(-1.0+z1);
		*(Gammaq.ptr+2) = *ptr21*(-1.0+z1);
		*(Gammaq.ptr+3) = z1*(*ptr12**ptr21*ps+*ptr22*(*ptr22+*ptr12**ptr21*ps-*ptr11**ptr22*ps)*pb);
	}
}

/* FT: Impuls- -> Ortsraum */
void FBT2 (t_fbt f) {
	unsigned i, m;

	double *ptr11_r, *ptr12_r, *ptr21_r, *ptr22_r;
	double *ptr11_q, *ptr12_q, *ptr21_q, *ptr22_q;

	for (i = 0, f.bessel.ptr = f.bessel.matrix, ptr11_r = f.mr.matrix, ptr12_r = ptr11_r+1, ptr21_r = ptr12_r+1, ptr22_r = ptr21_r+1;
			i < f.gitter; i++, ptr11_r+=4, ptr12_r+=4, ptr21_r+=4, ptr22_r+=4) {
		*ptr11_r = *ptr12_r = *ptr21_r = *ptr22_r = 0.0;
		for (m = 0, ptr11_q = f.mq.matrix, ptr12_q = ptr11_q+1, ptr21_q = ptr12_q+1, ptr22_q = ptr21_q+1;
				m < f.gitter; m++, f.bessel.ptr++, ptr11_q+=4, ptr12_q+=4, ptr21_q+=4, ptr22_q+=4) {
			*ptr11_r += *ptr11_q**f.bessel.ptr;
			*ptr12_r += *ptr12_q**f.bessel.ptr;
#ifndef SYMMETRIC_FT
			*ptr21_r += *ptr21_q**f.bessel.ptr;
#endif
			*ptr22_r += *ptr22_q**f.bessel.ptr;
		}
		*ptr11_r *= f.cons;
		*ptr12_r *= f.cons;
#ifdef SYMMETRIC_FT
		*ptr21_r = *ptr12_r;
#else
		*ptr21_r *= f.cons;
#endif
		*ptr22_r *= f.cons;
	}
	/* Fixme: Symmetrie ausnutzen */
}

/* Matrix symmetrisieren */
void symmetry (t_matrix matrix, t_matrix symmetric_matrix) {
	unsigned int i;

	matrix.ptr = matrix.matrix;
	symmetric_matrix.ptr = symmetric_matrix.matrix;
	for (i = 0; i < symmetric_matrix.size; i+=4, matrix.ptr+=4, symmetric_matrix.ptr+=4) {
		*(symmetric_matrix.ptr) = *matrix.ptr;
		*(symmetric_matrix.ptr+1) = 0.5*(*(matrix.ptr+1)+*(matrix.ptr+2));
		*(symmetric_matrix.ptr+2) = *(symmetric_matrix.ptr+1);
		*(symmetric_matrix.ptr+3) = *(matrix.ptr+3);
	}
}



/*******************************************
 * Routinen für die Zerlegungsbasiserzeugung
 *******************************************/
/* Multiindex-Matrix allokieren 1d */
double *alloc_matrix_1 (int anz_1) {   // Liefert einen Vektor mit Einträgen = 0
    	double *m;
    	unsigned int i;
    
    	m = (double *) malloc ((size_t) (anz_1 * sizeof (double)));
    	if (!m) 
		die ("allocation failure in alloc_matrix_1");

    	for (i = 0; i < anz_1; i++) 
	    	m[i] = 0.0;
    
    	return m;
}

/* Multiindex-Matrix allokieren 2d */
double **alloc_matrix_2(int anz_1,int anz_2) {
    	unsigned int i;
    	double **m;
    
    	m = (double **) malloc((size_t) (anz_1 * sizeof(double *)));
    	if (!m) 
		die ("allocation failure in alloc_matrix_2");
    	for (i = 0; i < anz_1; i++) 
		m[i] = alloc_matrix_1(anz_2);
   	return m;
}

/* Multiindex-Matrix freigeben 1d */
void free_matrix_1(double *m) { 
    	free(m);
}

/* Multiindex-Matrix freigeben 2d */
void free_matrix_2(double **m,int anz_1) { 
    	unsigned int i; 
    
    	for (i = 0; i < anz_1; i++)
		free_matrix_1(m[i]);
    	free(m);
}

/* Gegeben ein Vektor r_i und ein X, liefert dies i: |r_i-X|=min */
unsigned int find_i (double *r_, double X, unsigned int gitter) {
//#define __FIND_I__
	unsigned int i = gitter-1, j;
	double minimum = r_[i];
	
	for (j = 0; j < gitter; j++) {
#ifdef __FIND_I__
		out ("\ni:%d j:%d r_[j]:%lf X:%lf min:%lf, fabs:%lf", i, j, r_[j], X, minimum, fabs(r_[j]-X));
#endif
		if (fabs(r_[j]-X) <= minimum) {
			minimum = fabs(r_[j]-X);
			i = j;
		}
	}
#ifdef __FIND_I__
	out ("\n[find_i: r%Lf X%Lf min%Lf i%d]", r_[i], X, minimum, i);
#endif
	return i;
}

/* Invertieren einer 2x2-Matrix */
void invert (double **matrix, double **inverse, unsigned int dimension) {   
	int k, i, j;   // Counter
	double min_value = 0.0000000001;

	// Invertierbar hiermit?
	for (i = 0; i < dimension; i++) 
		if (matrix[i][i] < min_value) 
			die ("Invertierfehler: (Diagonalelemet %d) = %.64lf", i, matrix[i][i]);

	// Einheitsmatrix
	double **unity = alloc_matrix_2 (dimension, dimension);
	for (i = 0; i < dimension; i++) for (j = 0; j < dimension; j++) 
		unity[i][j] = (double)(!(i-j));

	// Invertierung
	double xmult;
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
	double sum;
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

/* Gegeben die Basispkt. r_i_a, Ortsmatrix r_, liefert es die Zerlegungsbasis P_, Q_ */
void pq_init (double r_i_a[], unsigned int nu, t_matrix r_, t_matrix P_, t_matrix Q_) {
//#define __PQ_INIT__
#ifdef __PQ_INIT__
	out ("\n=> PQ generieren");
#endif
	unsigned int i, m, n;
	unsigned int gitter = r_.size;

	// Stützpunkte wählen
	int i_[nu];   // Die nu i_alphas aus den Maximalstellen bestimmen
	for (n = 0; n < nu; n++) {
		i_[n] = find_i (r_.matrix, r_i_a[n], gitter);
#ifdef __PQ_INIT__
		out ("[%d %lf]", i_[n], r_i_a[n]);
#endif
			}
	
	// *** P_n^i konstruieren ***
#ifdef __PQ_INIT__
	out ("P_n^i ");
#endif
	double **P__ = alloc_matrix_2 (nu, gitter);
	for (n = 0; n < nu; n++)
		if (n == 0) {	// 1. Basis
			for (i = 0; i < gitter; i++) {
				P__[n][i] = 0.0;
				if (i <= i_[1]) P__[n][i] = (double)(i_[1]-i)/(double)(i_[1]);
			}
		} else if (n == 1) {	// 2. Basis
			for (i = 0; i < gitter; i++) {
				P__[n][i] = 0.0;
				if (0 <= i && i <= i_[1]) P__[n][i] = (double)(i)/(double)(i_[1]);
				else if (i_[1] <= i && i <= i_[2]) P__[n][i] = (double)(i_[2]-i)/(double)(i_[2]-i_[1]);
			}
		} else {	// Alle anderen Basen
			for (i = 0; i < gitter; i++) {
				P__[n][i] = 0.0;
				if (i_[n-1] <= i && i <= i_[n]) P__[n][i] = (double)(i-i_[n-1])/(-i_[n-1]+i_[n]);
				else if (i_[n] <= i && i <= i_[n+1]) P__[n][i] = (double)(i_[n+1]-i)/(i_[n+1]-i_[n]);
			}
		}

	// *** Q_n^i konstruieren ***
#ifdef __PQ_INIT__
	out ("Q_n^i ");
#endif
	double **Q__ = alloc_matrix_2 (nu, gitter);
	
	// Dazu die Matrix (PP)_mn = sum_i(P_m^i*P_n^i) berechnen
	double sum;
	double **PP = alloc_matrix_2 (nu, nu);
	for (m = 0; m < nu; m++) for (n = 0; n < nu; n++) {
		sum = 0.0;
		for (i = 0; i < gitter; i++)
			sum += P__[m][i]*P__[n][i];
		PP[m][n] = sum;
	}

	// Dazu die Inserve RR von PP berechnen
	double **RR = alloc_matrix_2 (nu, nu);
	invert (PP, RR, nu);

	// Und schließlich die Qs berechnen
	for (m = 0; m < nu; m++) 
		for (i = 0; i < gitter; i++) {
			Q__[m][i] = 0.0;
			for (n = 0; n < nu; n++) 
				Q__[m][i] += RR[m][n]*P__[n][i];
		}

	// P, Q zurückgeben in t_matrix
	for (i = 0, P_.ptr = P_.matrix, Q_.ptr = Q_.matrix; i < gitter; i++)
		for (n = 0; n < nu; n++, P_.ptr++, Q_.ptr++) {
			*P_.ptr = (double)P__[n][i];
			*Q_.ptr = (double)Q__[n][i];
#ifdef __PQ_INIT__
			out ("\n(p%lf q%lf)", P__[n][i], Q__[n][i]);
#endif
		}
			
}

/* Zerlegung einer Funktion in Grob-/Fein-Anteil */
void decomposite (t_composite d) {
//#define __DECOMPOSITE_DEBUG
	unsigned int i, m;

	// Grobanteil nullen
#ifdef __DECOMPOSITE_DEBUG
	printf ("\n=> Zerlegung: Grobanteil 1"); fflush (stdout);
#endif
	for (i = 0, d.coarse.ptr = d.coarse.matrix; i < d.coarse.size; i++, d.coarse.ptr++)
		*d.coarse.ptr = 0.0;

#ifdef __DECOMPOSITE_DEBUG
	printf (" 2"); fflush (stdout);
#endif
	// Grobanteil berechnen
	for (i = 0, d.function.ptr = d.function.matrix, d.Q_.ptr = d.Q_.matrix; i < d.gitter; i++, d.function.ptr+=4)
		for (m = 0, d.coarse.ptr = d.coarse.matrix; m < d.nu; m++, d.Q_.ptr++, d.coarse.ptr+=4) {
#ifdef __DECOMPOSITE_DEBUG
			printf (" (%f*%f)", *(d.function.ptr), *d.Q_.ptr); fflush (stdout);
#endif
			*(d.coarse.ptr) += *(d.function.ptr)**d.Q_.ptr;
			*(d.coarse.ptr+1) += *(d.function.ptr+1)**d.Q_.ptr;
			*(d.coarse.ptr+2) += *(d.function.ptr+2)**d.Q_.ptr;
			*(d.coarse.ptr+3) += *(d.function.ptr+3)**d.Q_.ptr;
		}

#ifdef __DECOMPOSITE_DEBUG
	printf (" Feinanteil 1"); fflush (stdout);
#endif
	// Feinanteil gleich Funktion...
	for (i = 0, d.fine.ptr = d.fine.matrix, d.function.ptr = d.function.matrix; i < d.fine.size; i++, d.fine.ptr++, d.function.ptr++)
		*d.fine.ptr = *d.function.ptr;

#ifdef __DECOMPOSITE_DEBUG
	printf (" 2"); fflush (stdout);
#endif
	// minus Summation: 
	for (i=0, d.fine.ptr=d.fine.matrix, d.P_.ptr=d.P_.matrix; i < d.gitter; i++, d.fine.ptr+=4)
		for (m = 0, d.coarse.ptr = d.coarse.matrix; m < d.nu; m++, d.coarse.ptr+=4, d.P_.ptr++) {
			*(d.fine.ptr) -= *(d.coarse.ptr)**d.P_.ptr;
			*(d.fine.ptr+1) -= *(d.coarse.ptr+1)**d.P_.ptr;
			*(d.fine.ptr+2) -= *(d.coarse.ptr+2)**d.P_.ptr;
			*(d.fine.ptr+3) -= *(d.coarse.ptr+3)**d.P_.ptr;
		}
}

/* Komposition von Grob-/Fein-Anteil zu einer Fkt. */
void composite (t_composite c) {
	/* Todo: Problem: 1. Pkt. immer == 0? Gilt auch beim monodispersen...*/
//#define _COMPOSITE_
	unsigned int i, m;

	// Zusammensetzen der Funktion
	for (i = 0, c.function.ptr = c.function.matrix, c.P_.ptr = c.P_.matrix, c.fine.ptr = c.fine.matrix; i < c.gitter; i++, c.fine.ptr+=4, c.function.ptr+=4){
#ifdef _COMPOSITE_
		if (i <= 3)
			out ("\ni:%d fine:%lf", i, *c.function.ptr);
#endif
		// Aus Fein-
		*(c.function.ptr) = *(c.fine.ptr);
		*(c.function.ptr+1) = *(c.fine.ptr+1);
		*(c.function.ptr+2) = *(c.fine.ptr+2);
		*(c.function.ptr+3) = *(c.fine.ptr+3);
		for (m = 0, c.coarse.ptr = c.coarse.matrix; m < c.nu; m++, c.coarse.ptr+=4, c.P_.ptr++) {
#ifdef _COMPOSITE_
			if (i <= 3)
				out (" cp:%lf", *c.coarse.ptr**c.P_.ptr);
#endif
			// und Grobanteil
			*(c.function.ptr) += *(c.coarse.ptr)**c.P_.ptr;
			*(c.function.ptr+1) += *(c.coarse.ptr+1)**c.P_.ptr;
			*(c.function.ptr+2) += *(c.coarse.ptr+2)**c.P_.ptr;
			*(c.function.ptr+3) += *(c.coarse.ptr+3)**c.P_.ptr;
		}
#ifdef _COMPOSITE_
		if (i <= 3)
			out (" -> %lf", *c.function.ptr);
#endif
	}
}

#ifdef _BASIS_
/* Debugausgabe der Basen */
void basis_out (t_basis_out b) {
	char filename[255];
	FILE *basis_file;
	unsigned int n, i;

	for (n = 0; n < b.nu; n++) {
		sprintf (filename, "P_%d.dat", n);
		basis_file = fopen (filename, "w");
		for (i = 0; i < b.gitter; i++) 
			fprintf (basis_file, "\n%f %f", b.r_.matrix[i], b.P_.matrix[i*b.nu+n]);
		fclose (basis_file);
	}
}
#endif


/***********
 * Marquardt
 ***********/
/* Ein PY/OZ-Zyklus im LM-Algorithmus*/
void lm_evaluate_default (double* coarse, int nu_, double* diff, void *data, int *info) {
	unsigned int i;
	lm_data_type *d = (lm_data_type*) data;

	// übergebener Grobanteil in comp, Zusammensetzen, PY, FBT, OZ, FBT, Zerlegung -> neuer Grobanteil -> Differenz
	for (i = 0, d->comp.coarse.ptr = d->comp.coarse.matrix; i < nu_; i++, d->comp.coarse.ptr++)
		*d->comp.coarse.ptr = coarse[i];
	
	composite (d->comp);

	// PY / HNC
#if APPROXIMATION == 1
	PY (d->Cr, d->Gammar, d->mayer);
#elif APPROXIMATION == 2
	HNC (d->Cr, d-> Gammar, d->Vr);
#endif
	
	// OZ
	d->fbt1.mr = d->Cr; d->fbt1.mq = d->Cq; 
	FBT1 (d->fbt1);	
	OZ (d->Cq, d->Gammaq, d->p_s, d->p_b);
	d->fbt2.mq = d->Gammaq; d->fbt2.mr = d->Gammar; 
	FBT2 (d->fbt2);
	decomposite (d->decomp);

	// Fixme: Symmetrisierung rein

	for (i = 0, d->comp.coarse.ptr = d->comp.coarse.matrix, d->decomp.coarse.ptr = d->decomp.coarse.matrix; i < nu_; i+=4, d->comp.coarse.ptr+=4, d->decomp.coarse.ptr+=4) {
		diff[i] = *d->comp.coarse.ptr-*d->decomp.coarse.ptr;
		diff[i+1] = *(d->comp.coarse.ptr+1)-*(d->decomp.coarse.ptr+1);
		diff[i+2] = *(d->comp.coarse.ptr+2)-*(d->decomp.coarse.ptr+2);
		diff[i+3] = *(d->comp.coarse.ptr+3)-*(d->decomp.coarse.ptr+3);
		if (isnan (diff[i]) || isnan (diff[i+1]) || isnan (diff[i+2]) || isnan (diff[i+3]))
			*info = -1;
	}
	
	*info = *info;	// -1 für Divergenz
}



/******************
 * Ausgabe-Routinen
 ******************/
/* Matrix in datei schreiben */
void write_matrix (t_matrix matrix, t_matrix vector, char *name, int pic) {
	char filename[255];
	sprintf (filename, "%s__pic%d.dat", name, pic);
	FILE *matrix_file = fopen (filename, "w");

	if (!matrix_file)
		die ("\nwrite_matrix: %s %d", name, pic);
	
	unsigned int i;
	
	for (i = 0, matrix.ptr = matrix.matrix, vector.ptr = vector.matrix; i < matrix.size; i+=4, matrix.ptr+=4, vector.ptr++) 
		fprintf (matrix_file, "%.16f %.16f %.16f %.16f %.16f\n", *vector.ptr, *(matrix.ptr), *(matrix.ptr+1), *(matrix.ptr+2), *(matrix.ptr+3));
	
	fclose (matrix_file);
}

/* Grob-/Feinanteil in Datei schreiben */
void coarse_fine_out (t_coarse_fine cf) {
	char filename_coarse[255], filename_fine[255];
	sprintf (filename_coarse, "coarse__pic%d", *cf.pic);
	sprintf (filename_fine, "fine__pic%d", *cf.pic);
	FILE *coarse_file = fopen (filename_coarse, "w"), *fine_file = fopen (filename_fine, "w");

	unsigned int i, m;
	
	for (i = 0; i < cf.gitter; i++) 
		fprintf (fine_file, "\n%.16f %.16f %.16f %.16f %.16f", cf.r_.matrix[i], 
				cf.fine.matrix[i*4+0*2+0], cf.fine.matrix[i*4+0*2+1], cf.fine.matrix[i*4+1*2+0], cf.fine.matrix[i*4+1*2+1]);
	
	for (m = 0; m < cf.nu; m++) 
		fprintf (coarse_file, "\n%.16f %.16f %.16f %.16f %.16f", cf.r_.matrix[i],
				cf.coarse.matrix[i*4+0*2+0], cf.coarse.matrix[i*4+0*2+1], cf.coarse.matrix[i*4+1*2+0], cf.coarse.matrix[i*4+1*2+1]);
	
	fclose (coarse_file);
	fclose (fine_file);
}

/* Alle Korrelationsmatrizen ausgeben */
void matrix_out (t_matrices m) {
	unsigned int i;
	t_matrix tmp = init_matrix (m.Gammar.size);

	// H(r) & H(q)
	for (i = 0, tmp.ptr = tmp.matrix, m.Gammar.ptr = m.Gammar.matrix, m.Cr.ptr = m.Cr.matrix; i < m.Cr.size; i++, m.Cr.ptr++, m.Gammar.ptr++, tmp.ptr++) 
		*tmp.ptr = *m.Gammar.ptr + *m.Cr.ptr;
	
	write_matrix (tmp, m.r_, "Hr", *m.pic);
	for (i = 0, tmp.ptr = tmp.matrix, m.Gammaq.ptr = m.Gammaq.matrix, m.Cq.ptr = m.Cq.matrix; i < m.Cq.size; i++, m.Cq.ptr++, m.Gammaq.ptr++, tmp.ptr++)
		*tmp.ptr = *m.Gammaq.ptr + *m.Cq.ptr;
	write_matrix (tmp, m.q_, "Hq", *m.pic);

	// Gamma(r) & Gamma(q)
	write_matrix (m.Gammar, m.r_, "Gammar", *m.pic);
	write_matrix (m.Gammaq, m.q_, "Gammaq", *m.pic);

	// C(r) & C(q)
	write_matrix (m.Cr, m.r_, "Cr", *m.pic);
	write_matrix (m.Cq, m.q_, "Cq", *m.pic);

	// S(q)
	for (i = 0, tmp.ptr = tmp.matrix, m.Gammaq.ptr = m.Gammaq.matrix, m.Cq.ptr = m.Cq.matrix; i < m.Cq.size; i+=4, m.Cq.ptr+=4, m.Gammaq.ptr+=4, tmp.ptr+=4) {
		*(tmp.ptr) = (*m.x_s) + (*m.x_s)*(*m.x_b)*(*m.p)*(*m.Gammaq.ptr+*m.Cq.ptr);
		*(tmp.ptr+1) = (*m.x_s)*(*m.x_b)*(*m.p)*(*(m.Gammaq.ptr+1)+*(m.Cq.ptr+1));
		*(tmp.ptr+2) = (*m.x_s)*(*m.x_b)*(*m.p)*(*(m.Gammaq.ptr+2)+*(m.Cq.ptr+2));
		*(tmp.ptr+3) = (*m.x_b) + (*m.x_s)*(*m.x_b)*(*m.p)*(*(m.Gammaq.ptr+3)+*(m.Cq.ptr+3));
	}
	write_matrix (tmp, m.q_, "Sq", *m.pic);

	// Aufräumen
	cleanup_matrix (tmp);
}



/***************************************
 * Interpolation mit bikubischen Splines
 ***************************************/
/* Kubischen Spline berechnen */
void spline (t_matrix x, t_matrix y, t_matrix y2) {
//#define __SPLINE__
	//double yp1 = 0.0, ypn = 0.0;	// 1st derivative at boundary points =0 => natural spline
	unsigned int i, k;
	double p, qn, sig, un;
#ifdef __SPLINE__
	printf ("=> Init"); fflush (stdout);
#endif
	t_matrix u = init_matrix (x.size-1);

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

/* Interpolation mit kubischem Spline */
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

/* Bikubischen Spline erstellen */
void splie2 (t_matrix x1a, t_matrix x2a, t_matrix ya, t_matrix y2a) {
	unsigned int j;
	t_matrix ya_tmp;
	t_matrix y2a_tmp;
	unsigned int m = x1a.size, n = x2a.size;
	ya_tmp.size = y2a_tmp.size = n;
	
//#define __SPLIE2__
#ifdef __SPLIE2__
	printf ("\n"); fflush (stdout);
#endif
	
	for (j = 0, ya.ptr = ya.matrix, y2a.ptr = y2a.matrix; j < m; j++, ya.ptr+=n, y2a.ptr+=n) {
		ya_tmp.matrix = ya.ptr;
		y2a_tmp.matrix = y2a.ptr;
#ifdef __SPLIE2__
		printf ("\nsplie2"); fflush (stdout);
#endif
		spline (x2a, ya_tmp, y2a_tmp);
#ifdef __SPLIE2__
		printf ("ok"); fflush (stdout);
#endif
	}
}

/* Interpolation mit bikubischem Spline */
double splin2 (t_matrix x1a, t_matrix x2a, t_matrix ya, t_matrix y2a, double x1, double x2) {
//#define __SPLIN2__
	double y;
	unsigned int j;

#ifdef __SPLIN2__
	printf ("\ninit"); fflush (stdout);
#endif
	unsigned int m = x1a.size;
	unsigned int n = x2a.size;
	t_matrix ytmp = init_matrix (m);	/* Fixme: Nur 1x */
	t_matrix yytmp = init_matrix (m);
	t_matrix ya_tmp, y2a_tmp;
	ya_tmp.size = y2a_tmp.size = n;

#ifdef __SPLIN2__
	printf (" for"); fflush (stdout);
#endif
	for (j = 0, ya.ptr = ya.matrix, y2a.ptr = y2a.matrix, yytmp.ptr = yytmp.matrix; j < m; j++, ya.ptr+=n, y2a.ptr+=n, yytmp.ptr++) {
#ifdef __SPLIN2__
		printf ("."); fflush (stdout);
#endif
		ya_tmp.matrix = ya.ptr;
		y2a_tmp.matrix = y2a.ptr;
		*yytmp.ptr = splint (x2a, ya_tmp, y2a_tmp, x2);
	}
#ifdef __SPLIN2__
	printf (" spline"); fflush (stdout);
#endif
	spline (x1a, yytmp, ytmp); //?
#ifdef __SPLIN2__
	printf (" splint"); fflush (stdout);
#endif
	y = splint (x1a, yytmp, ytmp, x1);

#ifdef __SPLIN2__
	printf (" cleanup"); fflush (stdout);
#endif
	cleanup_matrix (ytmp);
	cleanup_matrix (yytmp);
#ifdef __SPLIN2__
	printf (" ok"); fflush (stdout);
#endif
		
	return (y);
}



/*****************
 * Einleseroutinen
 *****************/
/* Matrix aus Datei lesen */
t_read_matrix read_matrix (t_matrix matrix, char *name) {
	FILE *matrix_file = fopen (name, "r");
	unsigned int i = 1;//, tmp;
	double rmax = 0;
	t_read_matrix result;

	matrix.ptr = matrix.matrix;
	//for (i = 0; i < matrix.size; i++, matrix.ptr++)
	//	fscanf (matrix_file, "%f %f %f %f %f\n", &rmax, (matrix.ptr), (matrix.ptr+1), (matrix.ptr+2), (matrix.ptr+3));

	result.stepsize = rmax/i;
	result.rows = i;
	fclose (matrix_file);

	return result;
}

/* Matrix interpolieren mit bikubischem Spline */
void interpol_matrix (t_matrix matrix, char *fileprefix, t_matrix q_, t_matrix n, double n_target) {
#define _INTERPOL_MATRIX_
	/* Cleanup: Variablennamen */
	// Übergeben wird:
	// - t_matrix n: Dichten zu denen S(q) vorliegen als: sq_0.n <- 6 stellen nach komma
	// - q_: gewünschte Impulsrasterung
	out ("\n-> S(q) extrapolieren");
	
	unsigned int i, j;
	unsigned int lines_in = 0;
	FILE *infile;
	char filename[255];

	// # Pkte. in einzulesenden Dateien bestimmen (!müssen homogen in allen Dateien sein!)
	out (" lines_in:");
	sprintf (filename, "%s_%lf", fileprefix, n.matrix[0]);
	infile = fopen (filename, "r");
	if (!infile)
		die ("startgamma: infile");
	double tmp;
	do {
		fscanf (infile, "\n%lf %lf %lf %lf %lf", &tmp, &tmp, &tmp, &tmp, &tmp);
		lines_in++;
	} while (!feof (infile));
	lines_in--;
	out ("%d", lines_in);
	fclose (infile);
	
	// Speicher allokieren
	out (" Init");
	t_matrix sq_in_11 = init_matrix (lines_in*n.size),
		 sq_in_12 = init_matrix (lines_in*n.size),
		 sq_in_21 = init_matrix (lines_in*n.size),
		 sq_in_22 = init_matrix (lines_in*n.size);
	t_matrix q_in = init_matrix (lines_in);
	t_matrix deriv_11 = init_matrix (lines_in*n.size),
		 deriv_12 = init_matrix (lines_in*n.size),
		 deriv_21 = init_matrix (lines_in*n.size),
		 deriv_22 = init_matrix (lines_in*n.size);
	t_matrix sqn_11 = init_matrix (q_.size),
		 sqn_12 = init_matrix (q_.size),
		 sqn_21 = init_matrix (q_.size),
		 sqn_22 = init_matrix (q_.size);
	
	// Bekannte S(q) einlesen
	out (" Matrix einlesen");
	for (i = 0, n.ptr = n.matrix; i < n.size; i++, n.ptr++) {
		out (" %s_%lf", fileprefix, *n.ptr);
		sprintf (filename, "%s_%lf", fileprefix, *n.ptr);
		infile = fopen (filename, "r");
		if (!infile)
			die ("\nstartgamma: infile\n");
		for (j = 0, sq_in_11.ptr = sq_in_11.matrix+i*lines_in, sq_in_12.ptr = sq_in_12.matrix+i*lines_in, 
				sq_in_21.ptr = sq_in_21.matrix+i*lines_in, sq_in_22.ptr = sq_in_22.matrix+i*lines_in, q_in.ptr = q_in.matrix; 
				j < lines_in; j++, sq_in_11.ptr++, sq_in_12.ptr++, sq_in_21.ptr++, sq_in_22.ptr++, q_in.ptr++)
			fscanf (infile, "\n%lf %lf %lf %lf %lf", q_in.ptr, sq_in_11.ptr, sq_in_12.ptr, sq_in_21.ptr, sq_in_22.ptr);
		fclose (infile);
	}

	// 2. Ableitung der eingelsenen S(q)
	out (" Ableitung");
	splie2 (n, q_in, sq_in_11, deriv_11);
	splie2 (n, q_in, sq_in_12, deriv_12);
	splie2 (n, q_in, sq_in_21, deriv_21);
	splie2 (n, q_in, sq_in_22, deriv_22);

	// Interpolieren 
	out (" Interpolieren: n=%lf", n_target);
	for (i = 0, sqn_11.ptr = sqn_11.matrix, sqn_12.ptr = sqn_12.matrix, sqn_21.ptr = sqn_21.matrix, sqn_22.ptr = sqn_22.matrix, q_.ptr = q_.matrix; 
			i < sqn_11.size; i++, sqn_11.ptr++, sqn_12.ptr++, sqn_21.ptr++, sqn_22.ptr++, q_.ptr++) 
		if (*q_.ptr < q_in.matrix[lines_in-1]) {
			*sqn_11.ptr = splin2 (n, q_in, sq_in_11, deriv_11, n_target, *q_.ptr);
			*sqn_12.ptr = splin2 (n, q_in, sq_in_12, deriv_12, n_target, *q_.ptr);
			*sqn_21.ptr = splin2 (n, q_in, sq_in_21, deriv_21, n_target, *q_.ptr);
			*sqn_22.ptr = splin2 (n, q_in, sq_in_22, deriv_22, n_target, *q_.ptr);
		} else  
			*sqn_11.ptr = *sqn_12.ptr = *sqn_21.ptr = *sqn_22.ptr = sqn_11.matrix[lines_in-2];

#ifdef _INTERPOL_MATRIX_
	out (" out");
	FILE *outfile = fopen ("interpol_matrix", "w");
	if (!outfile)
		die ("\nstartgamma: outfile\n");
	for (i = 0, q_.ptr = q_.matrix, sqn_11.ptr = sqn_11.matrix, sqn_12.ptr = sqn_12.matrix, sqn_21.ptr = sqn_21.matrix, sqn_22.ptr = sqn_22.matrix; 
			i < q_.size; i++, q_.ptr++, sqn_11.ptr++, sqn_12.ptr++, sqn_21.ptr++, sqn_22.ptr++) 
		fprintf (outfile, "\n%lf %lf %lf %lf %lf", *q_.ptr, *sqn_11.ptr, *sqn_12.ptr, *sqn_21.ptr, *sqn_22.ptr);
	fclose (outfile);
#endif

	// Cleanup
	out (" Cleanup");
	cleanup_matrix (sq_in_11);
	cleanup_matrix (sq_in_12);
	cleanup_matrix (sq_in_21);
	cleanup_matrix (sq_in_22);
	cleanup_matrix (q_in);
	cleanup_matrix (deriv_11);
	cleanup_matrix (deriv_12);
	cleanup_matrix (deriv_21);
	cleanup_matrix (deriv_22);
	cleanup_matrix (sqn_11);
	cleanup_matrix (sqn_12);
	cleanup_matrix (sqn_21);
	cleanup_matrix (sqn_22);
}

/* S(q) bikubisch interpolieren */
void interpol_sq (t_matrix Sq, t_matrix q_, t_matrix n, double n_target) {
#define _INTERPOL_SQ_
	// Übergeben wird:
	// - t_matrix n: Dichten zu denen S(q) vorliegen als: sq_0.n <- 6 stellen nach komma
	// - q_: gewünschte Impulsrasterung
	out ("\n-> S(q) extrapolieren");
	
	unsigned int i, j;
	unsigned int lines_in = 0;
	FILE *infile;
	char filename[255];

	// # Pkte. in einzulesenden Dateien bestimmen (!müssen homogen in allen Dateien sein!)
	out (" lines_in:");
	sprintf (filename, "sq_%lf", n.matrix[0]);
	infile = fopen (filename, "r");
	if (!infile)
		die ("startgamma: infile");
	double tmp;
	do {
		fscanf (infile, "\n%lf %lf %lf %lf %lf", &tmp, &tmp, &tmp, &tmp, &tmp);
		lines_in++;
	} while (!feof (infile));
	lines_in--;
	out ("%d", lines_in);
	fclose (infile);
	
	// Speicher allokieren
	out (" Init");
	t_matrix sq_in_11 = init_matrix (lines_in*n.size),
		 sq_in_12 = init_matrix (lines_in*n.size),
		 sq_in_21 = init_matrix (lines_in*n.size),
		 sq_in_22 = init_matrix (lines_in*n.size);
	t_matrix q_in = init_matrix (lines_in);
	t_matrix deriv_11 = init_matrix (lines_in*n.size),
		 deriv_12 = init_matrix (lines_in*n.size),
		 deriv_21 = init_matrix (lines_in*n.size),
		 deriv_22 = init_matrix (lines_in*n.size);
	t_matrix sqn_11 = init_matrix (q_.size),
		 sqn_12 = init_matrix (q_.size),
		 sqn_21 = init_matrix (q_.size),
		 sqn_22 = init_matrix (q_.size);
	
	// Bekannte S(q) einlesen
	out (" S(q) einlesen");
	for (i = 0, n.ptr = n.matrix; i < n.size; i++, n.ptr++) {
		out (" sq_%lf", *n.ptr);
		sprintf (filename, "sq_%lf", *n.ptr);
		infile = fopen (filename, "r");
		if (!infile)
			die ("\nstartgamma: infile\n");
		for (j = 0, sq_in_11.ptr = sq_in_11.matrix+i*lines_in, sq_in_12.ptr = sq_in_12.matrix+i*lines_in, 
				sq_in_21.ptr = sq_in_21.matrix+i*lines_in, sq_in_22.ptr = sq_in_22.matrix+i*lines_in, q_in.ptr = q_in.matrix; 
				j < lines_in; j++, sq_in_11.ptr++, sq_in_12.ptr++, sq_in_21.ptr++, sq_in_22.ptr++, q_in.ptr++)
			fscanf (infile, "\n%lf %lf %lf %lf %lf", q_in.ptr, sq_in_11.ptr, sq_in_12.ptr, sq_in_21.ptr, sq_in_22.ptr);
		fclose (infile);
	}

	// 2. Ableitung der eingelsenen S(q)
	out (" Ableitung");
	splie2 (n, q_in, sq_in_11, deriv_11);
	splie2 (n, q_in, sq_in_12, deriv_12);
	splie2 (n, q_in, sq_in_21, deriv_21);
	splie2 (n, q_in, sq_in_22, deriv_22);

	// Interpolieren 
	out (" Interpolieren: n=%lf", n_target);
	for (i = 0, sqn_11.ptr = sqn_11.matrix, sqn_12.ptr = sqn_12.matrix, sqn_21.ptr = sqn_21.matrix, sqn_22.ptr = sqn_22.matrix, q_.ptr = q_.matrix; 
			i < sqn_11.size; i++, sqn_11.ptr++, sqn_12.ptr++, sqn_21.ptr++, sqn_22.ptr++, q_.ptr++) 
		if (*q_.ptr < q_in.matrix[lines_in-1]) {
			*sqn_11.ptr = splin2 (n, q_in, sq_in_11, deriv_11, n_target, *q_.ptr);
			*sqn_12.ptr = splin2 (n, q_in, sq_in_12, deriv_12, n_target, *q_.ptr);
			*sqn_21.ptr = splin2 (n, q_in, sq_in_21, deriv_21, n_target, *q_.ptr);
			*sqn_22.ptr = splin2 (n, q_in, sq_in_22, deriv_22, n_target, *q_.ptr);
		} else  
			*sqn_11.ptr = *sqn_12.ptr = *sqn_21.ptr = *sqn_22.ptr = 1.0;

#ifdef _INTERPOL_SQ_
	out (" out");
	FILE *outfile = fopen ("interpol_sq", "w");
	if (!outfile)
		die ("\nstartgamma: outfile\n");
	for (i = 0, q_.ptr = q_.matrix, sqn_11.ptr = sqn_11.matrix, sqn_12.ptr = sqn_12.matrix, sqn_21.ptr = sqn_21.matrix, sqn_22.ptr = sqn_22.matrix; 
			i < q_.size; i++, q_.ptr++, sqn_11.ptr++, sqn_12.ptr++, sqn_21.ptr++, sqn_22.ptr++) 
		fprintf (outfile, "\n%lf %lf %lf %lf %lf", *q_.ptr, *sqn_11.ptr, *sqn_12.ptr, *sqn_21.ptr, *sqn_22.ptr);
	fclose (outfile);
#endif

	// Cleanup
	out (" Cleanup");
	cleanup_matrix (sq_in_11);
	cleanup_matrix (sq_in_12);
	cleanup_matrix (sq_in_21);
	cleanup_matrix (sq_in_22);
	cleanup_matrix (q_in);
	cleanup_matrix (deriv_11);
	cleanup_matrix (deriv_12);
	cleanup_matrix (deriv_21);
	cleanup_matrix (deriv_22);
	cleanup_matrix (sqn_11);
	cleanup_matrix (sqn_12);
	cleanup_matrix (sqn_21);
	cleanup_matrix (sqn_22);
}

/* C(q) aus S(q) berechnen */
void sq_to_cq (t_matrix sq, t_matrix cq, double n, double x_s, double x_b) {
	/* Todo: Geht noch nicht*/
	unsigned int i;

	double *sq11, *sq12, *sq21, *sq22;
	double *cq11, *cq12, *cq21, *cq22;
	double p1 = x_s*n, p2 = x_b*n,
	       x11 = x_s*x_s*n, x12 = x_s*x_b*n, x21 = x_b*x_s*n, x22 = x_b*x_b*n;

	for (i = 0, sq11 = sq.matrix, sq12 = sq11+1, sq21 = sq12+1, sq22 = sq21+1, cq11 = cq.matrix, cq12 = cq11+1, cq21 = cq12+1, cq22 = cq21+1;
			i < sq.size; i+=4, sq11+=4, sq12+=4, sq21+=4, sq22+=4, cq11+=4, cq12+=4, cq21+=4, cq22+=4) {
		*cq11 = (-p2*(-1.0+*sq11+*sq12**sq21)+p2*(-1.0+*sq11)**sq22-*sq21*x12+(-1.0+*sq11)*x22)/(p2*(-1.0+*sq22)*x11-p2**sq12*x21-x12*x21+x11*x22+
				p1*(-p2*(-1.0+*sq11+*sq12**sq21)+p2*(-1.0+*sq11)**sq22-*sq21*x12+(-1.0+*sq11)*x22));
		*cq12 = ((-1.0+*sq22)*x12-*sq12*x22)/(x12*x21+p2*(x11-*sq22*x11+*sq12*x21)-x11*x22+
				p1*(p2*(-1.0+*sq11+*sq12**sq21+*sq22-*sq11**sq22)+*sq21*x12+x22-*sq11*x22));
		*cq21 = (-*sq21*x11+(-1.0+*sq11)*x21)/(x12*x21+p2*(x11-*sq22*x11+*sq12*x21)-x11*x22+p1*(p2*(-1.0+*sq11+*sq12**sq21+*sq22-*sq11**sq22)+*sq21*x12+x22-*sq11*x22));
		*cq22 = (-p1*(-1.0+*sq11+*sq12**sq21)+p1*(-1.0+*sq11)**sq22+(-1.0+*sq22)*x11-*sq12*x21)/(p2*(-1.0+*sq22)*x11-p2**sq12*x21-x12*x21+x11*x22+p1*(-p2*(-1.0+*sq11
						+*sq12**sq21)+p2*(-1.0+*sq11)**sq22-*sq21*x12+(-1.0+*sq11)*x22));

	}
}

double delta_hs;
double delta_chi;

/* Mayerfunktion aus den Systemparametern berechnen */
void return_mayer (t_matrix mayer, t_matrix r_, t_matrix Vr, double density, double temperature, double b_field) {
	unsigned int i;
	
#if SYSTEM == 1	// Harte Scheibchen
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen
	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;
	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : 0.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : 0.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : 0.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : 0.0;
	}
	
	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : 0.0;
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : 0.0;
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : 0.0;
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : 0.0;
	}
#elif SYSTEM == 2	// Mix harter Scheibchen
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen
	double delta = delta_hs;
	r_b = 2.35 * pow (10.0, -6.0);
	r_s = r_b*delta;
	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;

	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : 0.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : 0.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : 0.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : 0.0;
	}
#elif SYSTEM == 3	// Dipolare harte Scheibchen
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen

	double kb = 1.3806505e-23;			// Boltzmannkonstante
	double mu_ = pow (10.0, -7.0);			// mu_0/4pi
	double chi_s = 6.6e-12, chi_b = 6.2e-11;	// magnetische Suszeptibilitäten der Teilchen
	double B = b_field;		 		// B-Feld
	double T = temperature;				// Temperatur	

#ifdef SYSTEM_3
	r_s = r_b = 1.4e-6;
	chi_s = chi_b = 6.2e-12;
#endif
	
	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;
	double norm = pow (sigma_bs, -3.0); //pow (sigma_bs, -3.0) * pow (density, 3.0/2.0);
		//pow (sigma_bs, 3.0) * pow (density, -3.0/2.0); //*8.0;
	// pow (sigma_bs, -3.0)

	double V_ss = (-mu_*B*B/(kb*T) * chi_s*chi_s * norm);
	double V_sb = (-mu_*B*B/(kb*T) * chi_s*chi_b * norm);
	double V_bs = (-mu_*B*B/(kb*T) * chi_b*chi_s * norm);
	double V_bb = (-mu_*B*B/(kb*T) * chi_b*chi_b * norm);

#define ATTRACT 0.001

	out ("\n-> reduziertes Potential: V_ss:%lf V_sb:%lf V_bs:%lf V_bb:%lf", V_ss, V_sb, V_bs, V_bb);

	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
#ifdef ATTRACT
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*pow(*r_.ptr,-3.0)-V_ss*pow(*r_.ptr,-3.0)*ATTRACT)-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*pow(*r_.ptr,-3.0)-V_bs*pow(*r_.ptr,-3.0)*ATTRACT)-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*pow(*r_.ptr,-3.0)-V_sb*pow(*r_.ptr,-3.0)*ATTRACT)-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*pow(*r_.ptr,-3.0)-V_bb*pow(*r_.ptr,-3.0)*ATTRACT)-1.0;
#else
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*pow(*r_.ptr,-3.0))-1.0;
#endif
	}

	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
#ifdef ATTRACT
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : (1.0-ATTRACT)*V_ss*pow(*r_.ptr,-3.0);
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*pow(*r_.ptr,-3.0);
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*pow(*r_.ptr,-3.0);
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*pow(*r_.ptr,-3.0);
#else
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : V_ss*pow(*r_.ptr,-3.0);
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*pow(*r_.ptr,-3.0);
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*pow(*r_.ptr,-3.0);
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*pow(*r_.ptr,-3.0);
#endif
	}
#elif SYSTEM == 4	// LJ
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen
	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;
	double kb = 1.3806505e-23;			// Boltzmannkonstante
	double mu_ = 10e-7;			// mu_0/4pi
	double T = temperature;

	double norm = 1.0;//pow (sigma_bs, -3.0);

	double V_ss = (-mu_/(kb*T) * norm);
	double V_sb = (-mu_/(kb*T) * norm);
	double V_bs = (-mu_/(kb*T) * norm);
	double V_bb = (-mu_/(kb*T) * norm);

	out ("\n-> reduziertes Potential: V_ss:%lf V_sb:%lf V_bs:%lf V_bb:%lf", V_ss, V_sb, V_bs, V_bb);

	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*(pow(*r_.ptr*sigma_bs,-6.0)-pow(*r_.ptr*sigma_bs,-12.0)))-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0)))-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0)))-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0)))-1.0;
	}
	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : V_ss*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0));
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0));
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0));
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*(pow(*r_.ptr*sigma_bs,-3.0)-pow(*r_.ptr*sigma_bs,-12.0));
	}
#elif SYSTEM == 5	// Exponentieller Abfall
	double V_ss = -50.0;//-120.0;
	double V_sb = -50.0;//30.0;
	double V_bs = -50.0;//;30.0;
	double V_bb = -50.0;//40.0;

	#define POTENTIAL ( exp (0.3/ *r_.ptr)-1.0 )
	out ("\n-> reduziertes Potential: V_ss:%lf V_sb:%lf V_bs:%lf V_bb:%lf", V_ss, V_sb, V_bs, V_bb);
	
	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.5) ? -1.0 : exp (V_ss*POTENTIAL)-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*POTENTIAL)-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*POTENTIAL)-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 0.5) ? -1.0 : exp (V_bb*POTENTIAL)-1.0;
	}

	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
		*(Vr.ptr) = (*(r_.ptr) < 1.5) ? -1.0/0.0 : V_ss*POTENTIAL;
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*POTENTIAL;
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*POTENTIAL;
		*(Vr.ptr+3) = (*(r_.ptr) < 0.5) ? -1.0/0.0 : V_bb*POTENTIAL;
	}
#undef POTENTIAL
#elif SYSTEM == 6
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen

	double kb = 1.3806505e-23;			// Boltzmannkonstante
	double mu_ = pow (10.0, -7.0);			// mu_0/4pi
	double chi_s = 6.6e-12, chi_b = 6.2e-11;	// magnetische Suszeptibilitäten der Teilchen
	double B = b_field;		 		// B-Feld
	double T = temperature;				// Temperatur	

#ifdef SYSTEM_3
	r_s = r_b = 1.4e-6;
	chi_s = chi_b = 6.2e-12;
#endif
	
	r_s = r_b*delta_hs;
	chi_s = chi_b*delta_chi;

	out (" MAYER: rs:%g rb:%g cs:%g cb:%g",r_s, r_b, chi_s, chi_b);
	
	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;
	double norm = pow (sigma_bs, -3.0); //pow (sigma_bs, -3.0) * pow (density, 3.0/2.0);
		//pow (sigma_bs, 3.0) * pow (density, -3.0/2.0); //*8.0;
	// pow (sigma_bs, -3.0)

	double V_ss = (-mu_*B*B/(kb*T) * chi_s*chi_s * norm);
	double V_sb = (-mu_*B*B/(kb*T) * chi_s*chi_b * norm);
	double V_bs = (-mu_*B*B/(kb*T) * chi_b*chi_s * norm);
	double V_bb = (-mu_*B*B/(kb*T) * chi_b*chi_b * norm);

//#define NEW_NORM
#ifdef NEW_NORM
	norm = pow (density/(sigma_bs*sigma_bs), -3.0/2.0);
	V_ss = (-mu_*B*B/(kb*T) * chi_s*chi_s);
	V_sb = (-mu_*B*B/(kb*T) * chi_s*chi_b);
	V_bs = (-mu_*B*B/(kb*T) * chi_b*chi_s);
	V_bb = (-mu_*B*B/(kb*T) * chi_b*chi_b);
#endif

	out ("\n-> reduziertes Potential: V_ss:%lf V_sb:%lf V_bs:%lf V_bb:%lf", V_ss, V_sb, V_bs, V_bb);

	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*pow(*r_.ptr,-3.0))-1.0;
#ifdef NEW_NORM
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*pow(*r_.ptr/norm,-3.0))-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*pow(*r_.ptr/norm,-3.0))-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*pow(*r_.ptr/norm,-3.0))-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*pow(*r_.ptr/norm,-3.0))-1.0;
#endif
	}

	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : V_ss*pow(*r_.ptr,-3.0);
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*pow(*r_.ptr,-3.0);
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*pow(*r_.ptr,-3.0);
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*pow(*r_.ptr,-3.0);
#ifdef NEW_NORM
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : V_ss*pow(*r_.ptr/norm,-3.0);
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*pow(*r_.ptr/norm,-3.0);
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*pow(*r_.ptr/norm,-3.0);
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*pow(*r_.ptr/norm,-3.0);
#endif
	}
#elif SYSTEM == 7
	double r_s = 1.4e-6, r_b = 2.35e-6;	// Radien der Teilchen
	double kb = 1.3806505e-23;			// Boltzmannkonstante
	double mu_ = pow (10.0, -7.0);			// mu_0/4pi
	double chi_s = 6.6e-12, chi_b = 6.2e-11;	// magnetische Suszeptibilitäten der Teilchen
	double B = b_field;		 		// B-Feld
	double T = temperature;				// Temperatur	

	double sigma_ss = r_s+r_s, sigma_bs = r_s+r_b, sigma_bb = r_b+r_b;
	double norm = pow (sigma_bs, -3.0); //pow (sigma_bs, -3.0) * pow (density, 3.0/2.0);

	double V_ss = (-mu_*B*B/(kb*T) * chi_s*chi_s * norm);
	double V_sb = (-mu_*B*B/(kb*T) * chi_s*chi_b * norm);
	double V_bs = (-mu_*B*B/(kb*T) * chi_b*chi_s * norm);
	double V_bb = (-mu_*B*B/(kb*T) * chi_b*chi_b * norm);


	out ("\n-> reduziertes Potential: V_ss:%lf V_sb:%lf V_bs:%lf V_bb:%lf", V_ss, V_sb, V_bs, V_bb);

	for (i = 0, mayer.ptr = mayer.matrix, r_.ptr = r_.matrix; i < mayer.size; i+=4, mayer.ptr+=4, r_.ptr++) {
		*(mayer.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0 : exp (V_ss*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_bs*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0 : exp (V_sb*pow(*r_.ptr,-3.0))-1.0;
		*(mayer.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0 : exp (V_bb*pow(*r_.ptr,-3.0))-1.0;
	}

	for (i = 0, Vr.ptr = Vr.matrix, r_.ptr = r_.matrix; i < Vr.size; i+=4, Vr.ptr+=4, r_.ptr++) {
		*(Vr.ptr) = (*(r_.ptr) < 1.0*sigma_ss/sigma_bs) ? -1.0/0.0 : V_ss*pow(*r_.ptr,-3.0);
		*(Vr.ptr+1) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_sb*pow(*r_.ptr,-3.0);
		*(Vr.ptr+2) = (*(r_.ptr) < 1.0) ? -1.0/0.0 : V_bs*pow(*r_.ptr,-3.0);
		*(Vr.ptr+3) = (*(r_.ptr) < 1.0*sigma_bb/sigma_bs) ? -1.0/0.0 : V_bb*pow(*r_.ptr,-3.0);
	}

#endif
}



/***************
 * Hauptschleife
 ***************/
int main (int argc, char* argv[]) {
	int a, b;//, c, d;   // Teilchen-(Matrix)-Counter
	int i, j;   // Orts-/Impuls-Counter
	int q;	// Impulscounter
	int m, n;   // Zerlegungsbasis-Counter / Impuls-Counter (FBT)


	/**********************
	 * Parameter einstellen
	 **********************/
	double p = 5e-2;;
	unsigned int gitter = 1*1024;
	double Rmax = 30.0;;
	double x_s = 0.3, x_b = 1.0-x_s;
#if SYSTEM == 1	
	/* Parameter für harte Scheibchen */
	out ("\n-> Harte Scheibchen");

#ifndef BASIS_LINE
	Rmax = 20.0;
	#define NU 20
	int nu = NU;
	double r_i_a[NU] = {.1, .2, .3, .4, .5, .65, .95, 1.2, 1.3, 1.5, 1.7, 2.0, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0, 5.0, 6.0};

	if (argc != 1+4)
		die ("\n<cmd> <dichte> <x_s> <rmax> <m>\n");
#else
	if ((int) argc <= 1+4)
		die ("\n<cmd> <dichte> <x_s> <rmax> <m> [basen...]\n");
	
#endif
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;
	Rmax = (double)(atof(argv[3]));
	gitter = (unsigned int)(atoi(argv[4]));

	double delta = 1.4/2.35;

	double temperature = 0.0;
	double b_field = 0.0;

#ifdef BASIS_LINE
	int NU = (int)argc - (1+4);
	int nu = NU;
	double *r_i_a = (double *) calloc ((size_t)nu,(size_t)(sizeof(double)));
	for (i = 0; i < NU; i++)
		r_i_a[i] = (double)(atof(argv[1+4+i]));
#endif

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf Rmax:%lf Pkte:%d", p, x_s, x_b, Rmax, gitter);
#elif SYSTEM == 2	// Mix harter Scheibchen
	/* Parameter für Mix harter Scheibchen */
	out ("\n-> Mix harter Scheibchen");

#ifndef BASIS_LINE
	Rmax = 30.0;
	#define NU 20
	int nu = NU;
	double r_i_a[NU] = {.1, .2, .3, .4, .5, .65, .95, 1.2, 1.3, 1.5, 1.7, 2.0, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0, 5.0, 6.0};

	if (argc != 1+5) 
		die ("\n<cmd> <dichte> <x1> <delta> <rmax> <m>\n");
#else
	if (argc <= 1+5)
		die ("\n<cmd> <dichte> <x1> <delta> <rmax> <m> [basen...]\n");

	int NU = (int)argc - (1+5);
	int nu = NU;
	double *r_i_a = (double *) calloc ((size_t)nu,(size_t)(sizeof(double)));
	for (i = 0; i < NU; i++)
		r_i_a[i] = (double)(atof(argv[1+5+i]));
#endif
	p = (double) atof (argv[1]);
	double xi = 0;
		//(double) atof (argv[2]);
	delta_hs =(double) atof (argv[3]);
	double delta = delta_hs;
	x_s = (double)atof(argv[2]);
		//xi /(1.0+xi);
	x_b = 1.0-x_s;
	gitter = (unsigned int)(atoi(argv[5]));
	Rmax = (double)(atof(argv[4]));

	double temperature = 0.0;
	double b_field = 0.0;
		
	out ("\n-> Parameter: p:%lf xi=x_s/x_b:%lf delta=r_s/r_b:%lf Rmax:%lf Pkte:%d", p, xi, delta, Rmax, gitter);
#elif SYSTEM == 3 	// Dipolare harte Scheibchen
	/* Parameter für dipolare harte Scheibchen */
	out ("\n-> Dipolare harte Scheibchen");
	
#ifndef BASIS_LINE
	Rmax = 100.0;	
	#define NU 45
	int nu = NU;
	double r_i_a[NU] = {.2, .65, .9, 1.4, 1.7, 2.0, 2.3, 2.7, 3.0, 3.3, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, // I1
	5.7, 6.0, 6.3, 6.6, 6.9, 7.2, 7.5,	// I2
	7.8, 8.1, 8.4, 8.7, 9.0, 9.5, 10.0, 10.5, 11.0, // I3
	11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.6, 16.0, // I4
	17.0, 18.0, // I5
	20.0};//, 36.0}; // I6
	//43.0, 49.0, 59.0, 70.0, 80.0, 90.0, 110.0, 130.0, 150.0, 175.0 }; // restliche
	//{.2, .5, .8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5.0, 5.3, 5.6, 5.9, 6.2, 6.5, 6.8, 7.1, 7.4, 7.7, 8.0, 8.3, 8.6, 8.9, 9.2, 9.5, 9.8, // 33
	//10.1, 10.4, 10.7, 11.0, 11.3, 11.6, 11.9, 12.2, 12.5, 12.8, 13.1, 13.4, 13.7, 14.0, 14.3, 14.6, 14.9, 15.5, 16.0, 16.0, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, //27
	//21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 35.0, 40.0}; //12
	//{.2, .65, .9, 1.2, 1.7, 2.3, 2.7, 3.3, 3.9, 4.6, 5.3, 6.0, 6.8, 7.7, 8.7, 9.8, 10.9, 11.9, 13.0, 14.0, 15.0, 16.5, 17.9, 19.8, 21.5, 23.0, 25.0, 27.5, 30.0, 33.0,
	//36.0, 40.0, 45.0, 51.0, 57.0};
	//double r_i_a[NU] = {.1, .25, .3, .4, .5, .65, .75, .9, 1.0, 1.2, 1.3, 1.5, 1.7, 2.0, 2.1, 2.3, 2.5, 2.7, 3.0, 4.5, 7.0, 9.0, 12.0, 15.0, 19.0, 22.0, 25.0, 28.0, 31.0, 40.0};
	//double r_i_a[NU] = {.1, .2, .45, .65, .9, 
	//	1.3, 1.9, 2.5, 3.0, 3.7, 
	//	4.0, 4.5, 5.2, 5.5, 6.0, 
	//	7.0, 8.5, 9.5, 10.0, 11.5, 
	//	12.5, 16.0, 18.0, 20.0, 25.0};

	gitter = .75*1024;
#endif
	
#ifndef BASIS_LINE
	if (argc != 1+6)
		die ("\n<cmd> <dichte [s^2]> <x_s> <temperatur [K]> <b [mT]> <rmax> <m>\n");
#else
	if ((int)argc <= 1+6)
		die ("\n<cmd> <dichte [s^2]> <x_s> <temperatur [K]> <b [mT]> <rmax> <m> [basen..]\n");
#endif
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;
	double temperature = (double) atof (argv[3]);
	double b_field = (double) (atof (argv[4])) * pow (10.0, -3.0);
	Rmax = (int)(atoi(argv[5]));
	gitter = (int)(atoi(argv[6]));

#ifdef BASIS_LINE
	int NU = (int)argc - (1+6);
	int nu = NU;
	double *r_i_a = (double *) calloc ((size_t)nu,(size_t)(sizeof(double)));
	for (i = 0; i < NU; i++)
		r_i_a[i] = (double)(atof(argv[1+6+i]));
#endif

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf T:%lf B:%lf Rmax:%lf Pkte:%d", p, x_s, x_b, temperature, b_field, Rmax, gitter);
#elif SYSTEM == 4	// LJ
	out ("\n-> Binäres LJ");

	Rmax = 20.0;
	#define NU 20
	int nu = NU;
	double r_i_a[NU] = {.1, .2, .3, .4, .5, .65, .95, 1.2, 1.3, 1.5, 1.7, 2.0, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0, 5.0, 6.0};
	gitter = 1*1024;

	if (argc != 1+3)
		die ("\n<cmd> <dichte> <x_s> <temperature>\n");
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;
	double temperature = (double) atof (argv[3]);

	double b_field = 0.0;

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf Rmax:%lf Pkte:%d", p, x_s, x_b, Rmax, gitter);
#elif SYSTEM == 5	// Exponentieller Abfall
	out ("\n-> Exponentieller Abfall");
	
	Rmax = 50.0;
	#define NU 20
	int nu = NU;
	double r_i_a[NU] = {.1, .2, .3, .4, .5, .65, .95, 1.2, 1.3, 1.5, 1.7, 2.0, 2.3, 2.5, 2.7, 3.0, 3.5, 4.0, 5.0, 6.0};
	gitter = 1024;

	if (argc != 1+2)
		die ("\n<cmd> <dichte> <x_s>\n");
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;

	double temperature = 0.0;
	double b_field = 0.0;

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf Rmax:%lf Pkte:%d", p, x_s, x_b, Rmax, gitter);
#elif SYSTEM == 6
	/* Parameter für dipolare harte Scheibchen */
	out ("\n-> Dipolare harte Scheibchen (mix)");
	
#ifndef BASIS_LINE
	Rmax = 100.0;	
	#define NU 45
	int nu = NU;
	double r_i_a[NU] = {.2, .65, .9, 1.4, 1.7, 2.0, 2.3, 2.7, 3.0, 3.3, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, // I1
	5.7, 6.0, 6.3, 6.6, 6.9, 7.2, 7.5,	// I2
	7.8, 8.1, 8.4, 8.7, 9.0, 9.5, 10.0, 10.5, 11.0, // I3
	11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.6, 16.0, // I4
	17.0, 18.0, // I5
	20.0};//, 36.0}; // I6
	//43.0, 49.0, 59.0, 70.0, 80.0, 90.0, 110.0, 130.0, 150.0, 175.0 }; // restliche
	//{.2, .5, .8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5.0, 5.3, 5.6, 5.9, 6.2, 6.5, 6.8, 7.1, 7.4, 7.7, 8.0, 8.3, 8.6, 8.9, 9.2, 9.5, 9.8, // 33
	//10.1, 10.4, 10.7, 11.0, 11.3, 11.6, 11.9, 12.2, 12.5, 12.8, 13.1, 13.4, 13.7, 14.0, 14.3, 14.6, 14.9, 15.5, 16.0, 16.0, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, //27
	//21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 35.0, 40.0}; //12
	//{.2, .65, .9, 1.2, 1.7, 2.3, 2.7, 3.3, 3.9, 4.6, 5.3, 6.0, 6.8, 7.7, 8.7, 9.8, 10.9, 11.9, 13.0, 14.0, 15.0, 16.5, 17.9, 19.8, 21.5, 23.0, 25.0, 27.5, 30.0, 33.0,
	//36.0, 40.0, 45.0, 51.0, 57.0};
	//double r_i_a[NU] = {.1, .25, .3, .4, .5, .65, .75, .9, 1.0, 1.2, 1.3, 1.5, 1.7, 2.0, 2.1, 2.3, 2.5, 2.7, 3.0, 4.5, 7.0, 9.0, 12.0, 15.0, 19.0, 22.0, 25.0, 28.0, 31.0, 40.0};
	//double r_i_a[NU] = {.1, .2, .45, .65, .9, 
	//	1.3, 1.9, 2.5, 3.0, 3.7, 
	//	4.0, 4.5, 5.2, 5.5, 6.0, 
	//	7.0, 8.5, 9.5, 10.0, 11.5, 
	//	12.5, 16.0, 18.0, 20.0, 25.0};

	gitter = .75*1024;
#elif SYSTEM == 7
	/* Parameter für dipolare harte Scheibchen */
	out ("\n-> Dipolare harte Scheibchen");
	
#ifndef BASIS_LINE
	Rmax = 100.0;	
	#define NU 45
	int nu = NU;
	double r_i_a[NU] = {.2, .65, .9, 1.4, 1.7, 2.0, 2.3, 2.7, 3.0, 3.3, 3.9, 4.2, 4.5, 4.8, 5.1, 5.4, // I1
	5.7, 6.0, 6.3, 6.6, 6.9, 7.2, 7.5,	// I2
	7.8, 8.1, 8.4, 8.7, 9.0, 9.5, 10.0, 10.5, 11.0, // I3
	11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.6, 16.0, // I4
	17.0, 18.0, // I5
	20.0};//, 36.0}; // I6
	//43.0, 49.0, 59.0, 70.0, 80.0, 90.0, 110.0, 130.0, 150.0, 175.0 }; // restliche
	//{.2, .5, .8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5.0, 5.3, 5.6, 5.9, 6.2, 6.5, 6.8, 7.1, 7.4, 7.7, 8.0, 8.3, 8.6, 8.9, 9.2, 9.5, 9.8, // 33
	//10.1, 10.4, 10.7, 11.0, 11.3, 11.6, 11.9, 12.2, 12.5, 12.8, 13.1, 13.4, 13.7, 14.0, 14.3, 14.6, 14.9, 15.5, 16.0, 16.0, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, //27
	//21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 35.0, 40.0}; //12
	//{.2, .65, .9, 1.2, 1.7, 2.3, 2.7, 3.3, 3.9, 4.6, 5.3, 6.0, 6.8, 7.7, 8.7, 9.8, 10.9, 11.9, 13.0, 14.0, 15.0, 16.5, 17.9, 19.8, 21.5, 23.0, 25.0, 27.5, 30.0, 33.0,
	//36.0, 40.0, 45.0, 51.0, 57.0};
	//double r_i_a[NU] = {.1, .25, .3, .4, .5, .65, .75, .9, 1.0, 1.2, 1.3, 1.5, 1.7, 2.0, 2.1, 2.3, 2.5, 2.7, 3.0, 4.5, 7.0, 9.0, 12.0, 15.0, 19.0, 22.0, 25.0, 28.0, 31.0, 40.0};
	//double r_i_a[NU] = {.1, .2, .45, .65, .9, 
	//	1.3, 1.9, 2.5, 3.0, 3.7, 
	//	4.0, 4.5, 5.2, 5.5, 6.0, 
	//	7.0, 8.5, 9.5, 10.0, 11.5, 
	//	12.5, 16.0, 18.0, 20.0, 25.0};

	gitter = .75*1024;
#endif
	
#ifndef BASIS_LINE
	if (argc != 1+6)
		die ("\n<cmd> <dichte [s^2]> <x_s> <temperatur [K]> <b [mT]> <rmax> <m>\n");
#else
	if ((int)argc <= 1+6)
		die ("\n<cmd> <dichte [s^2]> <x_s> <temperatur [K]> <b [mT]> <rmax> <m> [basen..]\n");
#endif
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;
	double temperature = (double) atof (argv[3]);
	double b_field = (double) (atof (argv[4])) * pow (10.0, -3.0);
	Rmax = (int)(atoi(argv[5]));
	gitter = (int)(atoi(argv[6]));

#ifdef BASIS_LINE
	int NU = (int)argc - (1+6);
	int nu = NU;
	double *r_i_a = (double *) calloc ((size_t)nu,(size_t)(sizeof(double)));
	for (i = 0; i < NU; i++)
		r_i_a[i] = (double)(atof(argv[1+6+i]));
#endif

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf T:%lf B:%lf Rmax:%lf Pkte:%d", p, x_s, x_b, temperature, b_field, Rmax, gitter);

#endif
	
#ifndef BASIS_LINE
	if (argc != 1+6)
		die ("\n<cmd> <dichte [s^2]> <x_s> <temperatur [K]> <b [mT]> <rmax> <m>\n");
#else
	if ((int)argc <= 1+8)
		die ("\n<cmd> <dichte [s^2]> <x_s> <delta_r> <delta_chi> <temperatur [K]> <b [mT]> <rmax> <m> [basen..]\n");
#endif
	p = (double) atof (argv[1]);
	x_s = (double) atof (argv[2]);
	x_b = 1.0-x_s;
	delta_hs = (double) atof(argv[3]);
	delta_chi = (double)atof(argv[4]);
	double delta = delta_hs;
	double temperature = (double) atof (argv[5]);
	double b_field = (double) (atof (argv[6])) * pow (10.0, -3.0);
	Rmax = (int)(atoi(argv[7]));
	gitter = (int)(atoi(argv[8]));

#ifdef BASIS_LINE
	int NU = (int)argc - (1+8);
	int nu = NU;
	double *r_i_a = (double *) calloc ((size_t)nu,(size_t)(sizeof(double)));
	for (i = 0; i < NU; i++)
		r_i_a[i] = (double)(atof(argv[1+8+i]));
#endif

	out ("\n-> Parameter: p:%lf x_s:%lf x_b:%lf T:%lf B:%lf Rmax:%lf Pkte:%d delta_r:%lf delta_chi:%lf", p, x_s, x_b, temperature, b_field, Rmax, gitter, delta_hs, delta_chi);

#endif


	//Packungsdichte!
#ifdef PACKING_FRACTION
	out ("\n-> Packingfraction: %lf", p);
	double delta = 0.59574468;
	double p_new = (p/(double)(M_PIl))*pow((delta+1.0),2.0)/(x_s*delta*delta+x_b);
	out (" -> ns_12^2=%lf", p_new);
	p = p_new;
#endif
	
	
	// Algorithmuseinstellungen ausgeben
#if APPROXIMATION == 1
	out ("\n-> PY");
#elif APPROXIMATION == 2
	out ("\n-> HNC");
#endif

#ifdef SYM
	printf ("\n-> Symmetriesiere Gamma(q) nach OZ"); fflush (stdout);
#if SYM == 1
	printf (" - alle 4 Matrixelemente gleich machen"); fflush (stdout);
#endif
#endif
#ifdef ASYMPTOTE
	printf ("\n-> Asymptote Gamma(r)=-1-C(r) r<sigma"); fflush (stdout);
#endif

	out ("\n-> %d Basen:", nu);
	for (i = 0; i < nu; i++)
		out (" %.2lf", r_i_a[i]);


	/************************
	 * Matrizen für Iteration
	 ************************/
	printf ("\n-> Matrizen initialisieren");
	t_matrix Cr = init_matrix (gitter*2*2);
	t_matrix Cq = init_matrix (gitter*2*2);
	t_matrix Gammar = init_matrix (gitter*2*2);
	t_matrix Gammaq = init_matrix (gitter*2*2);
#ifdef SYM
	t_matrix sym_Gammaq = init_matrix (gitter*2*2);
#endif
	t_matrix Vr = init_matrix (gitter*2*2);


	/******************
	 * Besselfunktionen
	 ******************/
	printf ("\n-> Besselfunktionen"); fflush (stdout);
	printf (" Fkt."); fflush (stdout);
	double Qmax = (double)(gsl_sf_bessel_zero_J0 (gitter+1)/Rmax);
	t_matrix r_ = init_matrix (gitter);
	for (i = 0, r_.ptr = r_.matrix; i < r_.size; i++, r_.ptr++)
		*r_.ptr = (double)(gsl_sf_bessel_zero_J0 (i+1)/Qmax);
	t_matrix q_ = init_matrix (gitter);
	for (m = 0, q_.ptr = q_.matrix; m < q_.size; m++, q_.ptr++)
		*q_.ptr = (double)(gsl_sf_bessel_zero_J0 (m+1)/Rmax);

	/* Fixme:
	 * Geht das in einer Struktur? */
	printf (" Fbt-Konstanten"); fflush (stdout);
	printf (" 1"); fflush (stdout);
	t_fbt fbt1 = {
		.cons = 4.0*(double)(M_PIl)/(Qmax*Qmax),
		.gitter = gitter,
		.bessel = init_matrix (gitter*gitter),
		.mr = Cr, 
		.mq = Cq
	};
	for (m = 0, fbt1.bessel.ptr = fbt1.bessel.matrix, q_.ptr = q_.matrix; m < q_.size; m++, q_.ptr++) 
		for (i = 0, r_.ptr = r_.matrix; i < r_.size; i++, fbt1.bessel.ptr++, r_.ptr++)
			*fbt1.bessel.ptr = j0(*q_.ptr**r_.ptr) / ( j1(*r_.ptr*Qmax)*j1(*r_.ptr*Qmax) );
	
	printf (" 2"); fflush (stdout);
	t_fbt fbt2 = {
		.cons = 1.0/((double)(M_PIl)*Rmax*Rmax),
		.gitter = gitter,
		.bessel = init_matrix (gitter*gitter),
		.mq = Gammaq,
		.mr = Gammar
	};
	for (i = 0, fbt2.bessel.ptr = fbt2.bessel.matrix, r_.ptr = r_.matrix; i < r_.size; i++, r_.ptr++) 
		for (m = 0, q_.ptr = q_.matrix; m < q_.size; m++, fbt2.bessel.ptr++, q_.ptr++)
			*fbt2.bessel.ptr = j0(*q_.ptr**r_.ptr) / ( j1(*q_.ptr*Rmax)*j1(*q_.ptr*Rmax) );

	
	/***************************
	 * Mayerfunktionen berechnen
	 ***************************/
	printf ("\n-> Mayerfunktion"); fflush (stdout);
	t_matrix mayer = init_matrix (gitter*2*2);
	return_mayer (mayer, r_, Vr, p, temperature, b_field);


	/************
	 * Startwerte
	 ************/
	printf ("\n-> Startwert: "); fflush (stdout);
#if STARTVALUE == 1	// Gamma(r) ist Startfkt.
	Gammar.ptr = Gammar.matrix;
	r_.ptr = r_.matrix;
	for (i = 0; i < Gammar.size; i+=4, Gammar.ptr+=4, r_.ptr++)
		*(Gammar.ptr) = *(Gammar.ptr+1) = *(Gammar.ptr+2) = *(Gammar.ptr+3) = 20.0*exp(-powl(*r_.ptr, 2.0)/4.0);
	printf (" Gamma(r) = 20*exp(-x^2/4)");
#elif STARTVALUE == 2   // Gamma(r) == 0
	for (i = 0, Gammar.ptr = Gammar.matrix; i < Gammar.size; i++, Gammar.ptr++)
		*Gammar.ptr = 0.0;
	// Ist nach Matrixinitialisierung immer der Fall
	printf (" Gamma(r) = 0");
#elif STARTVALUE == 3   // h(r) einlesen
	printf (" lese h(r) aus 'startvalue'"); fflush (stdout);
	t_read_matrix returnvalue;
	returnvalue = read_matrix (Gamma_q, gitter, "startvalue");
	printf (" %d Zeilen gelesen", returnvalue.rows);
	printf (" - konvertiere zu Gamma(r)"); fflush (stdout);
	for (i = 0; i < gitter; i++) for (a = 0; a < 2; a++) for (b = 0; b < 2; b++)
		Gamma_r[i][a][b] = Gamma_q[i][a][b];
#elif STARTVALUE == 4
	for (i = 0, Gammar.ptr = Gammar.matrix; i < Gammar.size; i++, Gammar.ptr++)
		*Gammar.ptr = 0.1;
#elif STARTVALUE == 5
	out (" Gamma(r) einlesen: ");
	read_matrix (Gammaq, "start_gamma");
	for (i = 0, Gammaq.ptr = Gammaq.matrix, Gammar.ptr = Gammar.matrix; i < Gammar.size; i++, Gammaq.ptr++, Gammar.ptr++)
		*Gammar.ptr = *Gammaq.ptr;
#endif	
	
	
	/*****************
	 * Zerlegungsbasis 
	 *****************/
	out ("\n-> Zerlegungsbasis");
	t_matrix P_ = init_matrix (gitter*nu),	// Zugriff P[i*nu+n]
		 Q_ = init_matrix (gitter*nu);	// Zugriff Q[i*nu+n]
	pq_init (r_i_a, nu, r_, P_, Q_);
//#define __BASIS_OUT__
#ifdef __BASIS_OUT__
	t_basis_out basis = {
		.P_ = P_,
		.gitter = gitter,
		.nu = nu,
		.r_ = r_
	};
	basis_out (basis);
#endif
//#define __TEST_PERP__
#ifdef __TEST_PERP__	// Auf Orthogonalität testn
	double sum;
	for (m = 0; m < nu; m++) {
		out ("\nm:%d   ", m);
		for (n = 0; n < nu; n++) {
			sum = 0.0;
			for (i = 0; i < gitter; i++) {
				sum += P_.matrix[i*nu+n]*Q_.matrix[i*nu+m];
			}
			out (" %lf", sum);
		}
	}
#endif

	// Strukturen für Zerlegung / Komposition erstellen
	printf (" Strukturen"); fflush (stdout);
	t_composite comp = {
		.function = Gammar,
		.coarse = init_matrix (nu*2*2),
		.fine = init_matrix (gitter*2*2),
		.P_ = P_,
		.Q_ = Q_,
		.gitter = gitter,
		.nu = nu
	};

	t_composite decomp = {
		.function = Gammar,
		.coarse = init_matrix (nu*2*2),
		.fine = init_matrix (gitter*2*2),
		.P_ = P_,
		.Q_ = Q_,
		.gitter = gitter,
		.nu = nu
	};


	/**************************
	 * Marquardt initialisieren
	 **************************/
#ifdef MARQUARDT
	out ("\n-> Marquardt initialisieren");
	
	lm_control_type control;
	lm_initialize_control (&control);
	
	lm_data_type data = {
		.Cr = Cr, .Cq = Cq, .Gammar = Gammar, .Gammaq = Gammaq,
		.p_s = x_s*p, .p_b = x_b*p,
		.mayer = mayer, .Vr = Vr,
		.fbt1 = fbt1, .fbt2 = fbt2,
		.comp = {
			.function = Gammar,
			.coarse = init_matrix (nu*2*2),
			.fine = init_matrix (gitter*2*2),
			.P_ = P_,
			.Q_ = Q_,
			.gitter = gitter,
			.nu = nu
		},
		.decomp = {
			.function = Gammar,
			.coarse = init_matrix (nu*2*2),
			.fine = init_matrix (gitter*2*2),
			.P_ = P_,
			.Q_ = Q_,
			.gitter = gitter,
			.nu = nu
		}
	};
	
#endif


	/***************
	 * Hauptschleife
	 ***************/
	out ("\n-> Hauptschleife");
	unsigned int pic = 0;
	unsigned int pic_max = 10000;
	double norm_fine = 0.0;  
	double norm_fine_max = 1e-12;

	time_t start, now;
	time (&start);

	out (" Ausgabe"); 
	// Ausgabe vorbereiten
	printf (" 1"); fflush (stdout);
	t_matrices matrices = {
		.Gammar = Gammar,
		.Gammaq = Gammaq,
		.Cr = Cr,
		.Cq = Cq,
		.r_ = r_,
		.q_ = q_,
		.pic = &pic,
		.x_s = &x_s,
		.x_b = &x_b,
		.p = &p
	};

	// Startwerte ausgeben
	printf (" 2"); fflush (stdout);
#ifdef _STARTVALUE_
	matrix_out (matrices);
#endif

	// Ausgabe von Grob- und Feinanteil vorbereiten
#ifdef _COARSE_FINE_
	printf (" 3"); fflush (stdout);
	t_coarse_fine coarse_fine = {
		.coarse = comp.coarse,
		.fine = comp.fine,
		.gitter = gitter,
		.nu = nu,
		.nr = &nr,
		.pic = &pic,
		.r_ = r_,
		.P_ = P_
	};
#endif


	/**************
	 * Testroutinen
	 **************/
#ifdef TEST
#if TEST == 1   // FBT 
	printf ("\n=> Teste FBT"); fflush (stdout);
	printf (" fbt1"); fflush (stdout);
	for (i = 0, Gammar.ptr = Gammar.matrix, r_.ptr = r_.matrix; i < gitter; i++, Gammar.ptr+=4, r_.ptr++) {
		*Gammar.ptr = exp (-*r_.ptr**r_.ptr);
		*(Gammar.ptr+1) = *(Gammar.ptr+2) = *(Gammar.ptr+3) = *(Gammar.ptr);
	}
	matrix_out (matrices);
	fbt1.mr = Gammar;
	fbt1.mq = Gammaq;
	FBT1 (fbt1);
	
	printf (" fbt2"); fflush (stdout);
	pic = 1;
	matrix_out (matrices);
	fbt2.mq = fbt1.mq;
	fbt2.mr = Cr;
	FBT2 (fbt2);

	printf (" 3"); fflush (stdout);
	pic = 2; 
	matrix_out (matrices);
#elif TEST == 2   // Zerlegung
	printf ("\n=> Teste Zerlegung");

	// Startwerte setzen & ausgeben
	printf (" Startwert"); fflush (stdout);
	pic = 1;
	for (i = 0, Gammar.ptr = Gammar.matrix, r_.ptr = r_.matrix; i < gitter; i++, Gammar.ptr+=4, r_.ptr++) 
		*Gammar.ptr = exp (-*r_.ptr**r_.ptr);
	matrix_out (matrices);

	// Zerlegen & ausgeben
	printf (" Zerlegen"); fflush (stdout);
	pic = 2;
	decomp.function = Gammar;
	decomposite (decomp);
	matrix_out (matrices);
	coarse_fine.coarse = decomp.coarse;
	coarse_fine.fine = decomp.fine;
	coarse_fine_out (coarse_fine);

	// Zusammensetzen & ausgeben
	printf (" Zusammensetzen"); fflush (stdout);
	pic = 3;
	for (i = 0, comp.coarse.ptr = comp.coarse.matrix, decomp.coarse.ptr = decomp.coarse.matrix; i < comp.coarse.size; i++, comp.coarse.ptr++, decomp.coarse.ptr++)
		*comp.coarse.ptr = *decomp.coarse.ptr;
	for (i = 0, comp.fine.ptr = comp.fine.matrix, decomp.fine.ptr = decomp.fine.matrix; i < comp.fine.size; i++, comp.fine.ptr++, decomp.fine.ptr++)
		*comp.fine.ptr = *decomp.fine.ptr;
	comp.function = Gammaq;
	composite (comp);
	matrix_out (matrices);
	coarse_fine.coarse = comp.coarse;
	coarse_fine.fine = comp.fine;
	coarse_fine_out (coarse_fine);
	printf (" ok"); fflush (stdout);
#elif TEST == 3	// Startgamma
	t_matrix Sq = init_matrix (gitter*4);
	t_matrix n_ = init_matrix (3);
	n_.matrix[0] = .95;
	n_.matrix[1] = .96;
	n_.matrix[2] = .97;
	double n_target = p = .965;
	
	interpol_sq (Sq, q_, n_, n_target);
	sq_to_cq (Sq, Cq, p, x_s, x_b);
	
	matrix_out (matrices);
	cleanup_matrix (n_);
	cleanup_matrix (Sq);
#elif TEST == 5	// Invertierung testen
	printf ("\n=> Invertierung"); fflush (stdout);

	// Matrix invertieren der Dimension dim
	unsigned int dim = 2;
	t_matrix unity = init_matrix (dim*dim);
	t_matrix inv = init_matrix (dim*dim);
	
	
	for (i = 0, unity.ptr = unity.matrix; i < unity.size; i++, unity.ptr++) {
		*unity.ptr = (double)(i+1);
		printf (" %f ", *unity.ptr);
	}
	
	printf ("\nMatrix vor Invertierung:");
	for (i = 0, unity.ptr = unity.matrix; i < dim; i++, unity.ptr++) {
		printf ("\n");
		for (j = 0; j < dim; j++) 
			printf ("%f ", unity.matrix[i*dim+j]);
	}

	/*
	for (i = 0, unity.ptr = unity.matrix; i < unity.size; i+=(dim+1), unity.ptr+=(dim+1))
		*unity.ptr = 2.0;
	*/

	/*
	unity.matrix[0*2+0]=1.0;
	unity.matrix[0*2+1]=2.0;
	unity.matrix[1*2+0]=3.0;
	unity.matrix[1*2+1]=4.0;
	*/
	
	invert (unity, inv, dim);

	printf ("\nInverse");
	for (i = 0, inv.ptr = inv.matrix; i < dim; i++, inv.ptr++) {
		printf ("\n");
		for (j = 0; j < dim; j++) 
			printf ("%.32f ", inv.matrix[i*dim+j]);
	}

	printf ("\nMatrix nach Invertierung:");
	for (i = 0, unity.ptr = unity.matrix; i < dim; i++, unity.ptr++) {
		printf ("\n");
		for (j = 0; j < dim; j++) 
			printf ("%f ", unity.matrix[i*dim+j]);
	}

	/* Test: Um absolut sicher zu sein: Bel. Matrix invertieren und dann mit der Inversen multiplizieren */

	cleanup_matrix (unity);
	cleanup_matrix (inv);
	printf ("ok"); fflush (stdout);
#elif TEST == 6
	out ("\n-> Zerlegung");
	
	for (i = 0, Gammar.ptr = Gammar.matrix, r_.ptr = r_.matrix; i < gitter; i++, Gammar.ptr+=4, r_.ptr++) 
		*Gammar.ptr = *(Gammar.ptr+1) = *(Gammar.ptr+2) = *(Gammar.ptr+3) = exp (-*r_.ptr**r_.ptr);
	
	matrix_out (matrices);
	decomposite (decomp);

	for (i = 0, decomp.coarse.ptr = decomp.coarse.matrix, comp.coarse.ptr = comp.coarse.matrix; i < decomp.coarse.size; i++, decomp.coarse.ptr++, comp.coarse.ptr++) {
		*comp.coarse.ptr = *decomp.coarse.ptr;
//		out ("\n%d coarse:%lf", i, *decomp.coarse.ptr);
	}

	for (i = 0, decomp.fine.ptr = decomp.fine.matrix, comp.fine.ptr = comp.fine.matrix; i < decomp.fine.size; i++, decomp.fine.ptr++, comp.fine.ptr++) {
		*comp.fine.ptr = *decomp.fine.ptr;
//		out ("\n%d fine:%lf", i, *decomp.fine.ptr);
	}

	i = 0;
	out ("\n--> %d %lf", i, Gammar.matrix[i]);
	
	pic++;
	composite (comp);
	out ("\n--> %d %lf", i, Gammar.matrix[i]);
	matrix_out (matrices);

	out ("\n--> %d %ld", i, Gammar.matrix[i]);
#elif TEST == 7	// Startgamma einlesen
	t_matrix n_ = init_matrix (3);

	// Fkt. einlesen bei den Dichten:
	n_.matrix[0] = .95;
	n_.matrix[1] = .96;
	n_.matrix[2] = .97;
	double n_target = p = .965;	
	
	out ("\n=> Startgamma einlesen");
	char fileprefix[255];
	sprintf (fileprefix, "gammar");
	interpol_matrix (Gammar, fileprefix, r_, n_, n_target);	

	matrix_out (matrices);
	cleanup_matrix (n_);
#endif
	// Ende der Testroutinen
	goto clean;
#endif


	double out_period_norm = 1e-1;
	unsigned int out_period = 1;
	
	/*******************************
	 * Grobanteil dynamisch anpassen
	 *******************************/
#ifdef DYN_COARSE
	out ("\n-> Grobanteilgenauigkeit adaptiv anpassen");
	// Schrittweitenparameter
	double dyn_coarse_step = 1e-1,
		dyn_coarse = 1e-3,
		dyn_coarse_max = 1e-12;
	double norm_fine_max_dyn_coarse = norm_fine_max*10;
	// Startwerte
	control.epsilon = dyn_coarse;
	control.ftol = dyn_coarse;
	control.xtol = dyn_coarse;
	control.gtol = dyn_coarse;
	norm_fine_max = dyn_coarse_max;
#else
	control.epsilon = control.ftol = control.xtol = control.gtol = 1e-15;
#endif

#ifdef FINE_RESCALE
	double alpha = 1.0;
#endif
	
	// Erste Zerlegung
	printf ("\n-> Erste Zerlegung:"); fflush (stdout);
	decomposite (decomp); 
	for (i = 0, decomp.fine.ptr = decomp.fine.matrix, comp.fine.ptr = comp.fine.matrix; i < decomp.fine.size; i++, decomp.fine.ptr++, comp.fine.ptr++)
		*comp.fine.ptr = *decomp.fine.ptr;
	for (i = 0, decomp.coarse.ptr = decomp.coarse.matrix, comp.coarse.ptr = comp.coarse.matrix; i < decomp.coarse.size; i++, decomp.coarse.ptr++, comp.coarse.ptr++)
		*comp.coarse.ptr = *decomp.coarse.ptr;

#ifdef T_PATH	// Temperaturwanderung
	double t_step = 0.001;
	double t_max = 10000.0, t_min = temperature;
	temperature = t_max;
	do {
	out ("\n-> *** T=%lf ***", temperature);
	return_mayer (mayer, r_, Vr, p, temperature, b_field);
#endif
	double norm_fine_old = norm_fine;
	
	/***********
	 * Marquardt
	 ***********/
	do {
		pic++;
		out ("\n%d ", pic);

		// Grob
		for (i = 0, data.comp.fine.ptr = data.comp.fine.matrix, comp.fine.ptr = comp.fine.matrix; i < comp.fine.size; i++, data.comp.fine.ptr++, comp.fine.ptr++)	
			// Feinanteil für LM fest vorgeben; Grobanteil wird ja eben variiert
			*data.comp.fine.ptr = *comp.fine.ptr;
		lm_minimize (4*nu, 4*nu, comp.coarse.matrix, lm_evaluate_default, lm_print_default, &data, &control);	/* Todo: Stimmt die Übergabe? */

		// Fein
		norm_fine_old = norm_fine;
		norm_fine = 0.0;
		for (i = 0, comp.fine.ptr = comp.fine.matrix, data.decomp.fine.ptr = data.decomp.fine.matrix; i < comp.fine.size; i++, comp.fine.ptr++, data.decomp.fine.ptr++) {
			norm_fine += pow ((*comp.fine.ptr-*data.decomp.fine.ptr), 2.0);//(*comp.fine.ptr-*data.decomp.fine.ptr)*(*comp.fine.ptr-*data.decomp.fine.ptr);
			*comp.fine.ptr = *data.decomp.fine.ptr;
		}
		
		if (norm_fine != 0) 
			norm_fine = sqrt (norm_fine);
/*
#ifdef FINE_RESCALE
		if (norm_fine > 1.0)
			alpha = 1.0/norm_fine;
		else
			alpha = 0.9;

		for (i = 0, comp.fine.ptr = comp.fine.matrix, data.decomp.fine.ptr = data.decomp.fine.matrix; i < comp.fine.size; i++, comp.fine.ptr++, data.decomp.fine.ptr++) {
			*comp.fine.ptr = alpha**data.decomp.fine.ptr+(1.0-alpha)**comp.fine.ptr;
		}
#endif
*/		
		/*
		if (norm_fine_old < norm_fine) {
			out ("!!");
			for (i = 0; i < gitter; i++) 
				comp.fine.matrix[i] = .1*data.decomp.fine.matrix[i] + (1.0-.1)*comp.fine.matrix[i];
		} else {
			for (i = 0; i < gitter; i++)
				comp.fine.matrix[i] = data.decomp.fine.matrix[i];
		}
		*/
#ifdef DEBUG
		matrix_out (matrices);
#endif

		/* Ausgabe */
		if (norm_fine < out_period_norm) {
			matrix_out (matrices);
			out_period_norm = .1*norm_fine;
			out_period *= 10;
		} 
		if ((pic % out_period) == 0) 
			matrix_out (matrices);

		out (" fine: %g", norm_fine);

#ifdef DYN_COARSE
		if (norm_fine < dyn_coarse && dyn_coarse > dyn_coarse_max) {//norm_fine_max_dyn_coarse && dyn_coarse > dyn_coarse_max) {
			dyn_coarse *= dyn_coarse_step;
			control.epsilon = dyn_coarse;
			control.ftol = dyn_coarse;
			control.xtol = dyn_coarse;
			control.gtol = dyn_coarse;
			out ("\n-> Grobanteilgenauigkeit anpassen %g\n", dyn_coarse);
		}
#endif
		time (&now);
		out (" time: %.0lf", difftime (now, start));
#ifdef DYN_COARSE
	} while (norm_fine > norm_fine_max && pic < pic_max);
#else
	} while (norm_fine > norm_fine_max && pic < pic_max);
#endif

	
#ifdef T_PATH	// Temperaturwanderung
	matrix_out (matrices);
	temperature -= t_step;
	} while (temperature > t_min);
#endif
			

	a = b = j = q = n = 1;	// Gegen Warnungen
#ifdef DYN_COARSE
	norm_fine_max_dyn_coarse = 0.0;
#endif
	
	/***********
	 * Aufräumen
	 ***********/
	goto clean;	// Gegen unused label
clean:
	out ("\n-> Aufräumen"); 
	matrix_out (matrices);
	
	out (" Korrelatoren");
	cleanup_matrix (Cr); out (".");
	cleanup_matrix (Cq); out (".");
#ifdef SYM
	cleanup_matrix (sym_Gammaq);
#endif
	cleanup_matrix (mayer); out (".");
	cleanup_matrix (Gammar); out (".");
	cleanup_matrix (Gammaq); out (".");

	cleanup_matrix (Vr);
	
	out (" Bessel");
	cleanup_matrix (r_);
	cleanup_matrix (q_);
	cleanup_matrix (fbt1.bessel);
	cleanup_matrix (fbt2.bessel);
	
	out (" Basis");
	cleanup_matrix (P_);
	cleanup_matrix (Q_);
	
	out (" Grob/Fein");
	cleanup_matrix (comp.fine);
	cleanup_matrix (comp.coarse);
	cleanup_matrix (decomp.fine);
	cleanup_matrix (decomp.coarse);

	out (" Marquardt");
	cleanup_matrix (data.comp.fine);
	cleanup_matrix (data.comp.coarse);
	cleanup_matrix (data.decomp.fine);
	cleanup_matrix (data.decomp.coarse);

#ifdef BASIS_LINE
	free (r_i_a);
#endif
	
	
	printf ("\n-> Fertig\n");
	return 0;
}
