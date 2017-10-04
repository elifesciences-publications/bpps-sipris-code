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

#if !defined(_SQS_TYP_)
#define _SQS_TYP_

#include "stdinc.h"
#include "alphabet.h"
#include "residues.h"
#include "sequence.h"

//******************** Add these later to sequence.h *******************
e_type	*MkEconoSeq(Int4 TotalSeqs, unsigned char **S,char **info, unsigned short *lengths);
void	NilEconoSeq(e_type);
//******************** Add these later to sequence.h *******************


class sqs_typ{
public:
		sqs_typ(){ assert(!"Illegal constructor"); }
  		sqs_typ(FILE *,a_type); // read in binary sqs_typ from file
  		sqs_typ(char *,a_type); // open file and read fasta sequences...
  		~sqs_typ(){ Free(); }
  void			Put(FILE *fp, double cutoff);
  void			Write(FILE *fp, double cutoff);
  e_type        Seq(Int4 i) {
			return MkSeq(seq_info[i],seq_len[i],seq_ptr[i]);
                }
  unsigned char *PtrSeq(UInt4 i){
                   if(i <= num_seqs) return seq_ptr[i]; else return 0;
                }
  Int4          SeqId(UInt4 i){
                   if(i <= num_seqs) return seq_id[i]; else return 0;
                }
  Int4 		NmbrSqs() { return num_seqs; }
private:
  Int4          MaxLength;
  Int4          MinLength;
  unsigned char	*HugeSeqMem;	// Huge memory array for seq_ptr;
  char		*HugeInfoMem;	// Huge memory array for seq_info;
  Int4		num_seqs;
  unsigned char **seq_ptr;	// pointer to sequences.
  char		**seq_info;	// pointer to sequence deflines.
  unsigned short *seq_len;
  e_type	*E;
  Int4		total_residues;
  Int4		*counts;
  Int4		*seq_id;
  char		**seq_start;
  char		**seq_ident;
  void		Init(char *,Int4);
  void          init_memory(Int4);
  void		Free();
  void		ReadFasta();
  void		WriteFasta();
  // void		ReadBinary();	// Wait on this for now.
  // void		WriteBinary();	// Wait on this for now.
  a_type	AB;
 };
#endif

