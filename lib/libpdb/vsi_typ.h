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

/* vsi_typ.h - CHAIN analysis 'View Structural Interaction' program. */
#if !defined (VSI_TYP)
#define VSI_TYP
#include <time.h>
#include "stdinc.h"
#include "afnio.h"
#include "ras_typ.h"
#include "pml_typ.h"
#include "mheap.h"


//**************************** VSI Type ****************************

struct trace_type {	// 27-55A.Y100
	UInt4 start;		// 27
	UInt4 end;		// 55
        char    chain;          // chain
        char    backbone;          // color 
        char    color;          // color 
        Int4	width;          // wireframe width
};

struct molecule_type {	// ATG802A_n1-h.X100
        char    *mol;        // molecule (e.g., 'ATG')
        Int4    id;          // molecule id (e.g., 802.
        char    chain;       // chain
        char    *atom;     // atom type (e.g., 'od1')
        char    *hydrogen; // attached hydrogen (e.g., 'od1')
        char    color;       // color 
        char    clouds;          // dot cloud
        Int4	width;          // wireframe width
};

struct vsi_residue_type {	// R127A_nh2-2hh2.X
        char    res;            // residue type (e.g., 'T')
        Int4    site;           // residue site.
        char    chain;          // chain
        char    *atom;        // atom type (e.g., 'od1')
        char    *hydrogen;    // attached hydrogen (e.g., 'od1')
        char    color;          // color 
        char    clouds;          // dot cloud
        Int4	width;          // wireframe width
};

union item_type {
	vsi_residue_type	residue;
        molecule_type		molecule;
        trace_type 		trace;
        Int4			*view;	// view: file18:R=-109,-45,109;T=-20,-12,506.
        char    		*cmd;	// echo; clear; file; File
};

// commands: echo clear File file 

struct item_list_type {
	item_type	item;
	char	 	type;
	Int4		file;
	char		*comment;
	item_list_type	*next;
};

typedef item_list_type *item_list_typ;

typedef item_list_type *vsi_typ;

struct dots_list_type {
	char		*str;	// e.g., A138D (was defined as "ala138D.CB")
	char	  	color;	// color {X} --> Ion Dots; else sidechain dots
	dots_list_type	*next;
};

dots_list_type *MakeDotsList(char *str,char color);

// Constructors:
item_list_type *MakeItemNodeTrace(UInt4 start, UInt4 end,
                char chain, char backbone, char color,Int4 width);
item_list_type *MakeItemNodeMol(char *mol, Int4 id, char chain, char *atom, 
		char *hydrogen, char color, char clouds, Int4 width);
item_list_type *MakeItemNodeRes(char res, Int4 site, char chain, char *atom,
		char *hydrogen, char color, char clouds, Int4 width);
item_list_type *MakeItemNode(char *cmd,char type);
item_list_type *MakeItemNodeView(Int4 *view);

// change settings:

// Information:
Int4    NumberAtomsItem(item_type item,char type);

//************************* VSI Heap Type ************************
typedef struct {
        mh_type mH;             /** minmax heap for best hits **/
        item_list_typ *VSI;    /** vsi item **/
} vsiheap_type;

typedef vsiheap_type *vsih_typ;
// typedef vsiheap_type *vsh_typ;
/*********************************************************************/

/******************************* private *******************************/

/******************************* PUBLIC *******************************/
vsih_typ  MakeVsiHeap(Int4 hpsz);
void    NilVsiHeap(vsih_typ H);
Int4    InsertVsiHeap(item_list_typ VSI, double key, vsih_typ H);
item_list_typ DelMinVsiHeap(double *key, Int4 *Item, vsih_typ H);
item_list_typ DelMaxVsiHeap(double *key, Int4 *Item, vsih_typ H);

/*********************************************************************/
#define nVsiHeap(H)     (ItemsInMheap((H)->mH))

//=============== vsi_typ.cc ======================
Int4    PutCmdItem(FILE *fp,item_type item,char type,Int4 file);
Int4    PutViewItem(FILE *fp,item_type item,char type,Int4 key_file);
Int4    PutMolItem(FILE *fp,item_type item,Int4 file);
Int4    PutTraceItem(FILE *fp,item_type item,Int4 file);
Int4    PutResItem(FILE *fp,item_type item,Int4 add,Int4 file);
Int4    PrintItems(FILE *fp, item_list_type *list,char *filename, char *FileComment,
		char key_chain, Int4 key_file,Int4 add,Int4 file,BooLean verbose);
Int4    PrintItems(FILE *fp, item_list_type *list,Int4 n, char **filename, 
		char **FileComment,char *key_chain, Int4 *key_file,Int4 add,Int4 file,
		BooLean verbose);
item_list_type *SortItems(item_list_type *list);
item_list_type *SortItems(item_list_type *list,BooLean purge);
item_list_typ SortItemsByChains(item_list_typ list);

void    FreeNode(item_list_typ node);
BooLean IdenticalNodes(item_list_typ node1, item_list_typ node2);
dots_list_type *MakeDotsList(char *str,char color);

//=============== vsi_ras.cc ======================
Int4	PrintItemsRasmol(FILE *fp,BooLean look_in_datadir, item_list_type *list,
		char *filename, char key_chain, BooLean cis_trans_pro,
		Int4 wire_width,Int4 spacefill);
Int4    PrintDotsItems(FILE *fp,item_list_type *list,ras_typ *ras);

//=============== vsi_pml.cc ======================
Int4    DefineItemsPyMol(FILE *fp, item_list_type *list, char key_chain, pml_typ *pml);
Int4	PrintItemsPyMol(FILE *fp,BooLean look_in_datadir, item_list_type *list,
		char *filename, char key_chain, BooLean cis_trans_pro,
		Int4 wire_width,Int4 spacefill);
Int4    PrintDotsItems(FILE *fp,item_list_type *list,pml_typ *pml);

#endif

