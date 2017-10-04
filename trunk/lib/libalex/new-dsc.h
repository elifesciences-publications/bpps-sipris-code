/******************************************************************************************
    Copyright (C) 1997-2014 Andrew F. Neuwald, Cold Spring Harbor Laboratory
    and the University of Maryland School of Medicine.

    Permission is hereby granted, free of charge, to any person obtaining a copy of 
    this software and associated documentation files (the "Software"), to deal in the 
    Software without restriction, including without limitation the rights to use, copy, 
    modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the 
    following conditions:

    The above copyright notice and this permission notice shall be included in all 
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
    PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
    OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
    OTHER DEALINGS IN THE SOFTWARE.

    For further information contact:
         Andrew F. Neuwald
         Institute for Genome Sciences and
         Department of Biochemistry & Molecular Biology
         University of Maryland School of Medicine
         801 West Baltimore St.
         BioPark II, Room 617
         Baltimore, MD 21201
         Tel: 410-706-6724; Fax: 410-706-1482; E-mail: aneuwald@som.umaryland.edu
 ******************************************************************************************/

#if !defined(_NEW_DSC_)
#define _NEW_DSC_

#include "cmsa.h"

# define Max_hom_seqs 5000
# define Max_length 5000
# define Max_median 10
# define LINE_LENGTH 5000
# define MAX_NAME_LENGTH 4500

typedef struct {
  	double 	infoa;
	double 	infob;
	double 	infoc;
	double 	edge_dist;
	double 	deletion;
	double 	insertion;
	double 	hydro_a;
	double 	hydro_b;
	double 	cons_a;
	double 	cons_b;
	double 	s_infoa;
	double 	s_infob;
	double 	s_infoc;
	double 	s_edge_dist;
	double 	s_deletion;
	double 	s_insertion;
	double 	s_hydro_a;
	double 	s_hydro_b;
	double 	s_cons_a;
	double 	s_cons_b;
	double 	prob_a;		/* probability of alpha-helical secondary structure */
	double 	prob_b;		/* probability of beta-strand secondary structure */
	double 	prob_c;		/* probability of random coil secondary structure */
	char 	prediction;
} dsc_att_vector;

typedef struct {
	int	error_line;
	int	sequence_length;
	char	**names;
	char 	**sequence;
	int 	hom_len;
	int 	filter_level;	
	int 	rem_isolated;
	int 	clean;
	int 	clean_length;
	int 	clean_percent;
	int 	s;
	dsc_att_vector	*att_array;
	Int4	NumHomSeq; // AFN addition...
} dsc_type;
typedef dsc_type *dsc_typ;

/******************************** PUBLIC ***********************************/
double **ComputeDSC(FILE *outf, double *predicted_accuracy, FILE *fp);
double **ComputeDSC(FILE *outf, double *predicted_accuracy, Int4 block, 
	Int4 left_flank, Int4 right_flank, cma_typ cma);
/******************************** PRIVATE **********************************/
dsc_typ	MkDSC(FILE *fp);
dsc_typ MkDSC(Int4 block, unsigned char left_flank, unsigned char right_flank,cma_typ cma);
void	NilDSC(dsc_typ D);
int PredictSeqDSC(dsc_typ D);
int PredictGorDSC(int pos1, double *infa, double *infb, double *infc, dsc_typ D);
int estimate_probs_DSC(dsc_typ D);
double pred_acc_DSC(dsc_typ D);
int form_res_atts_DSC(dsc_typ D);
int discrim1_DSC(double *ratio_a, double *ratio_b, dsc_typ D);
int discrim2_DSC(double ratio_a, double ratio_b, 
	double res_ratio_l, double res_ratio_w, 
	double res_ratio_a, double res_ratio_g, 
	double res_ratio_y, dsc_typ D);
int smooth_all_DSC(dsc_typ D);
int filter_DSC(dsc_typ D);
int remove_isolated_DSC(dsc_typ D);
double pcd(double tot_ian, double tot_ibn, double tot_icn);        
double res_ratio(char residue_type, dsc_typ D);
int copy_info(int flag, int info_type,
double  work_array[Max_length], dsc_typ D);
int median3(double work_array1[Max_length+1], double work_array2[Max_length+1]);
int median5(double work_array1[Max_length+1], double work_array2[Max_length+1]);
double m3(int pos, double work_array[Max_length+1]);
int median2(double work_array1[Max_length+1], double work_array2[Max_length+1]);
double m2(int pos, double work_array[Max_length+1]);
int median4(double work_array1[Max_length+1], double work_array2[Max_length+1]);
int clear_w(double w[Max_median]);
int insert(double p1, double w[Max_median]);
double medi(int m_length, double w[Max_median]);
int endpoints(double work_array[Max_length]);                        
int hanning(double work_array1[Max_length + 1], double work_array2[Max_length + 1]);
int rough(double work_array1[Max_length + 1], double work_array2[Max_length + 1], 
	double work_array3[Max_length + 1]);                        
int add(double work_array1[Max_length + 1], double work_array2[Max_length + 1], 
	double work_array3[Max_length + 1]);
double calc_edge_dist(int i, dsc_typ D);
double predict_del(int pos1, dsc_typ D);   
double predict_ins(int pos1, dsc_typ D);
double predict_h_moment(int pos1, double alpha, dsc_typ D);
double predict_c_moment(int pos, double rot_angle, dsc_typ D);
int get_res_eis(int hom, int pos, double n, double a,
	double *sin2, double *cos, dsc_typ D);
int get_res_con(int pos, double n, double a,double *sin2, double *cos2, dsc_typ D);
char get_res1(int hom, int n1,dsc_typ D);
double get_cons(int n1, dsc_typ D);  
double get_eisenberg(char nc0);
double entrop(int res_no, int no_seqs);
char pc(double tot_ian, double tot_ibn, double tot_icn);
double discrima1(int i, dsc_typ D);
double discrimb1(int i, dsc_typ D);
double discrimc1(int i, dsc_typ D);
double discrima2(int i, double ratio_a, double ratio_b,
double res_ratio_l, double res_ratio_w, double res_ratio_a, double res_ratio_g, 
	double res_ratio_y, dsc_typ D); 
double discrimb2(int i, double ratio_a, double ratio_b, double res_ratio_l, double res_ratio_w, 
	double res_ratio_a, double res_ratio_g, double res_ratio_y, dsc_typ D);
double discrimc2(int i, double ratio_a, double ratio_b, double res_ratio_l, 
	double res_ratio_w, double res_ratio_a, double res_ratio_g, double res_ratio_y, dsc_typ D);
int get_discrete(double diff);
int f1a(int i, dsc_typ D);
int f1b(int i, dsc_typ D);
int f2a(int i, dsc_typ D);        
int f3a(int i, dsc_typ D);
int f4a(int i, dsc_typ D);
int f5a(int i, dsc_typ D);
int f6a(int i, dsc_typ D);
int f7a(int i, dsc_typ D);
int f8a(int i, dsc_typ D);
int f9a(int i, dsc_typ D);
int f10a(int i, dsc_typ D);
int neighbour(int i, char *cn1, char *cn2, char *cn3, char *cp1, char *cp2,
	char *cp3, dsc_typ D);
int acces(char c, int pos, int *ian, int *ibn, int *icn);
int pos_to_pos(int pos1);
int aa_no(char c);
int edit_data(int read_length, dsc_typ D);
int upper_lower(int read_length, dsc_typ D);
int edit_probe_inserts(int read_length, dsc_typ D);
int calc_insert_lengths(int insert_position, dsc_typ D);
int edit_seqs(int insert_position, int insert_length, dsc_typ D);
int edit_begin(int read_length, dsc_typ D);
int remove_stars(dsc_typ D);
int e_beginnings(dsc_typ D);
int e_ends(dsc_typ D);
int get_seq_length(dsc_typ D);
int clean_up(int seq_length, dsc_typ D);
int dump(int seq_length, dsc_typ D);
int read_clustalw_sequence(FILE *fp, dsc_typ D);
int find_clustalw_start(FILE *fp_in, FILE *fp_out, dsc_typ D);
int read_clustalw_block(FILE *fp_in, FILE *fp_out, char first[LINE_LENGTH],
        int block, int *out_sequence_length, int *file_end, dsc_typ D);
int read_clustalw_line(char temp[LINE_LENGTH],char name[MAX_NAME_LENGTH],
        int block, int homologous_sequences, int out_sequence_length, dsc_typ D);
int read_clustalw_seq(char temp[LINE_LENGTH], int seq_no,
        int out_sequence_length, int i, dsc_typ D);

#endif

