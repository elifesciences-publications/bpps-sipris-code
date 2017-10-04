/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
#line 34 "chn_vsi.y"

#include "chn_vsi.h"
/**************************** Global Variables ******************************/
int	yyparse();
int     yylex();
char	*FileName=0;
char	*PDB_File=0;
Int4	NUMBER_OF_PDB_FILES;
char	*PDB_File_Comment=0;
Int4	KEY_FILE=0;
char	KEY_CHAIN=0;
UInt4	PARSELINE=1;

Int4    TRACE_WIDTH=60;
Int4    WIRE_WIDTH=60;
Int4    THIN_WIRE_WIDTH=50;
Int4    SPACEFILL=100;
Int4    BIG_SPACEFILL=160;
BooLean	CIS_TRANS_PRO=FALSE;
BooLean	USE_CURRENT_DIR=FALSE;
BooLean	VERBOSE=FALSE;
BooLean	debug=FALSE;

extern int yyerror(char *s);

/****************************************************************************/
vsi_typ	HEAD_NODE=0;
vsi_typ	TAIL_NODE=0;
/********************* Declarations of Action Symbols ***********************/

/* Line 371 of yacc.c  */
#line 99 "chn_vsi.tab.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "chn_vsi.tab.h".  */
#ifndef YY_YY_CHN_VSI_TAB_H_INCLUDED
# define YY_YY_CHN_VSI_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     RESIDUE = 258,
     MOLECULE = 259,
     FILENAME = 260,
     PDBFILE = 261,
     LETTER = 262,
     COLOR = 263,
     CHAIN = 264,
     ECHO_CMD = 265,
     COMMENT = 266,
     AN_ATOM = 267,
     HYDROGEN = 268,
     INTEGER = 269,
     VIEW = 270,
     CARBONYL = 271,
     JUNK = 272,
     NEWLINE = 273,
     END_OF_FILE = 274,
     CLEAR = 275
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 83 "chn_vsi.y"

		int	dummy;		/* no (null) attributes */
		char 	*pcval;		/* pointer to character string */
		char 	*mval;		/* pointer to molecule character string */
		char 	cval;		/* character */
		Int4	*Vval;		/* view of protein */
		char 	aval[9];	/* atom type*/
		Int4	ival;		/* integer type attribute */
		double	rval;		/* double (real) type attribute */
		vsi_typ	Lval;	// item list 
		char	Cval[3];	// color [0]='Y'; [1]=dots? (boolean)
		char	Sval[200];	// color [0]='Y'; [1]=dots? (boolean)
	

/* Line 387 of yacc.c  */
#line 177 "chn_vsi.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_CHN_VSI_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 205 "chn_vsi.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  15
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   143

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  34
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  64
/* YYNRULES -- Number of states.  */
#define YYNSTATES  131

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   275

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      32,    33,     2,    25,    21,    22,    24,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    27,     2,
       2,     2,     2,     2,    28,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    23,    29,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    30,     2,    31,    26,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     7,    12,    15,    20,    24,    27,    31,
      34,    36,    40,    44,    49,    52,    57,    62,    68,    70,
      73,    75,    78,    80,    82,    85,    87,    89,    91,    93,
      95,    97,   103,   110,   112,   114,   116,   118,   120,   122,
     127,   133,   140,   148,   155,   163,   172,   182,   187,   193,
     200,   208,   217,   227,   229,   231,   235,   239,   241,   244,
     246,   248,   251,   255,   258
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      35,     0,    -1,    52,    37,    41,    -1,    11,    52,    37,
      41,    -1,    37,    41,    -1,    11,    52,    36,    52,    -1,
      52,    36,    52,    -1,    36,    52,    -1,     1,    52,    41,
      -1,    37,    50,    -1,    50,    -1,    37,    38,    52,    -1,
      37,    38,    21,    -1,    37,    38,    11,    52,    -1,    51,
      52,    -1,    51,    52,    38,    52,    -1,    51,    52,    38,
      21,    -1,    51,    52,    38,    11,    52,    -1,    45,    -1,
      45,    14,    -1,    46,    -1,    46,    14,    -1,    39,    -1,
      43,    -1,    43,    14,    -1,    42,    -1,    20,    -1,    40,
      -1,    10,    -1,    19,    -1,    15,    -1,    14,    22,    14,
      44,     8,    -1,    14,    22,    14,     9,    44,     8,    -1,
      23,    -1,    24,    -1,    25,    -1,    26,    -1,    27,    -1,
      28,    -1,     3,    14,    24,    49,    -1,     3,    14,     9,
      24,    49,    -1,     3,    14,    29,    47,    24,    49,    -1,
       3,    14,     9,    29,    47,    24,    49,    -1,     3,    14,
      29,    16,    24,    49,    -1,     3,    14,     9,    29,    16,
      24,    49,    -1,     3,    14,    29,    47,    22,    48,    24,
      49,    -1,     3,    14,     9,    29,    47,    22,    48,    24,
      49,    -1,     4,    14,    24,    49,    -1,     4,    14,     9,
      24,    49,    -1,     4,    14,    29,    47,    24,    49,    -1,
       4,    14,     9,    29,    47,    24,    49,    -1,     4,    14,
      29,    47,    22,    48,    24,    49,    -1,     4,    14,     9,
      29,    47,    22,    48,    24,    49,    -1,    12,    -1,    13,
      -1,    30,     8,    31,    -1,    32,     8,    33,    -1,     8,
      -1,    50,    17,    -1,    17,    -1,     6,    -1,     6,    11,
      -1,    52,    11,    18,    -1,    52,    18,    -1,    18,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,   103,   103,   106,   109,   112,   113,   114,   115,   118,
     119,   122,   125,   128,   131,   132,   133,   134,   137,   138,
     139,   140,   141,   142,   143,   144,   147,   148,   150,   152,
     154,   156,   158,   162,   163,   164,   165,   166,   167,   169,
     171,   173,   175,   178,   181,   185,   188,   193,   195,   197,
     199,   201,   203,   208,   210,   212,   213,   214,   217,   219,
     221,   222,   225,   226,   227
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "RESIDUE", "MOLECULE", "FILENAME",
  "PDBFILE", "LETTER", "COLOR", "CHAIN", "ECHO_CMD", "COMMENT", "AN_ATOM",
  "HYDROGEN", "INTEGER", "VIEW", "CARBONYL", "JUNK", "NEWLINE",
  "END_OF_FILE", "CLEAR", "','", "'-'", "'^'", "'.'", "'+'", "'~'", "':'",
  "'@'", "'_'", "'{'", "'}'", "'('", "')'", "$accept", "start",
  "fatal_error", "item_list", "item", "cmd", "echo", "end_file", "view",
  "trace", "backbone", "residue", "molecule", "atom", "hydrogen", "color",
  "junk", "pdbfile", "nl", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,    44,    45,    94,    46,    43,   126,    58,    64,    95,
     123,   125,    40,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    34,    35,    35,    35,    35,    35,    35,    35,    36,
      36,    37,    37,    37,    37,    37,    37,    37,    38,    38,
      38,    38,    38,    38,    38,    38,    39,    39,    40,    41,
      42,    43,    43,    44,    44,    44,    44,    44,    44,    45,
      45,    45,    45,    45,    45,    45,    45,    46,    46,    46,
      46,    46,    46,    47,    48,    49,    49,    49,    50,    50,
      51,    51,    52,    52,    52
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     3,     4,     2,     4,     3,     2,     3,     2,
       1,     3,     3,     4,     2,     4,     4,     5,     1,     2,
       1,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     5,     6,     1,     1,     1,     1,     1,     1,     4,
       5,     6,     7,     6,     7,     8,     9,     4,     5,     6,
       7,     8,     9,     1,     1,     3,     3,     1,     2,     1,
       1,     2,     3,     2,     1
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    60,     0,    59,    64,     0,     0,     0,    10,
       0,     0,     0,    61,     0,     1,     7,     0,     0,    28,
       0,    30,    29,    26,     0,    22,    27,     4,    25,    23,
      18,    20,     9,    58,    14,     0,    63,     0,     0,     8,
       0,     0,     0,     0,     0,     0,    12,    11,    24,    19,
      21,     0,    62,     6,     2,     5,     3,     0,     0,     0,
       0,     0,     0,     0,    13,     0,    16,    15,     0,     0,
      57,     0,     0,    39,    53,     0,     0,     0,     0,    47,
       0,     0,    33,    34,    35,    36,    37,    38,     0,    17,
      40,     0,     0,     0,     0,     0,     0,     0,    48,     0,
       0,     0,     0,    31,     0,     0,     0,    55,    56,    43,
      54,     0,    41,     0,     0,     0,    49,    32,    44,     0,
      42,     0,     0,    50,     0,     0,    45,     0,    51,    46,
      52
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     6,     7,     8,    24,    25,    26,    27,    28,    29,
      88,    30,    31,    76,   111,    73,     9,    10,    11
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -95
static const yytype_int8 yypact[] =
{
      83,   -15,    30,   -15,   -95,   -95,    57,   -15,    76,    70,
     -15,    31,    43,   -95,    31,   -95,    28,    74,    78,   -95,
      81,   -95,   -95,   -95,     2,   -95,   -95,   -95,   -95,    93,
      96,   108,    70,   -95,    94,   106,   -95,   -15,    76,   -95,
     -15,    76,     0,     9,   111,   -15,   -95,    28,   -95,   -95,
     -95,     4,   -95,    28,   -95,    28,   -95,    26,    -4,    53,
      27,    -4,   114,    50,    28,   -15,   -95,    28,    -4,    55,
     -95,    98,   119,   -95,   -95,   104,    46,    -4,   114,   -95,
      59,    92,   -95,   -95,   -95,   -95,   -95,   -95,   121,    28,
     -95,   107,    89,   101,    97,    -4,   120,    -4,   -95,    99,
     120,    -4,   126,   -95,    -4,   120,    -4,   -95,   -95,   -95,
     -95,   112,   -95,   120,    -4,   113,   -95,   -95,   -95,   115,
     -95,    -4,   116,   -95,    -4,    -4,   -95,    -4,   -95,   -95,
     -95
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -95,   -95,    71,    88,   109,   -95,   -95,   -11,   -95,   -95,
      54,   -95,   -95,   -57,   -94,   -61,    -6,   -95,     7
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      79,    39,    32,     5,    70,    80,   115,    90,    12,    57,
      14,   119,    92,    45,    16,    65,    98,    34,    60,   122,
       5,    99,     5,    46,    58,    66,    71,    54,    72,    59,
      56,    47,    32,    61,   109,    32,   112,     2,    62,    35,
     116,    13,    35,   118,    53,   120,    36,    55,     4,    36,
      68,    77,    64,   123,    35,    69,    78,    15,    67,    81,
     126,    36,    22,   128,   129,    74,   130,    74,    96,    75,
      97,    91,    89,    82,    83,    84,    85,    86,    87,    17,
      18,   100,    37,   101,     1,    40,    19,    33,    42,     2,
      20,    21,    43,     4,     3,    22,    23,    17,    18,    38,
       4,     5,    41,    44,    19,    35,    93,    48,    20,    21,
      49,   105,    36,   106,    23,    82,    83,    84,    85,    86,
      87,   113,    50,   114,    52,    63,    74,    94,    95,   103,
     108,   104,   107,   110,   117,   102,   121,   124,     0,   125,
     127,     0,     0,    51
};

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-95)))

#define yytable_value_is_error(Yytable_value) \
  YYID (0)

static const yytype_int8 yycheck[] =
{
      61,    12,     8,    18,     8,    62,   100,    68,     1,     9,
       3,   105,    69,    11,     7,    11,    77,    10,     9,   113,
      18,    78,    18,    21,    24,    21,    30,    38,    32,    29,
      41,    24,    38,    24,    95,    41,    97,     6,    29,    11,
     101,    11,    11,   104,    37,   106,    18,    40,    17,    18,
      24,    24,    45,   114,    11,    29,    29,     0,    51,     9,
     121,    18,    19,   124,   125,    12,   127,    12,    22,    16,
      24,    16,    65,    23,    24,    25,    26,    27,    28,     3,
       4,    22,    11,    24,     1,    14,    10,    17,    14,     6,
      14,    15,    14,    17,    11,    19,    20,     3,     4,    11,
      17,    18,    14,    22,    10,    11,     8,    14,    14,    15,
      14,    22,    18,    24,    20,    23,    24,    25,    26,    27,
      28,    22,    14,    24,    18,    14,    12,     8,    24,     8,
      33,    24,    31,    13,     8,    81,    24,    24,    -1,    24,
      24,    -1,    -1,    34
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     6,    11,    17,    18,    35,    36,    37,    50,
      51,    52,    52,    11,    52,     0,    52,     3,     4,    10,
      14,    15,    19,    20,    38,    39,    40,    41,    42,    43,
      45,    46,    50,    17,    52,    11,    18,    36,    37,    41,
      36,    37,    14,    14,    22,    11,    21,    52,    14,    14,
      14,    38,    18,    52,    41,    52,    41,     9,    24,    29,
       9,    24,    29,    14,    52,    11,    21,    52,    24,    29,
       8,    30,    32,    49,    12,    16,    47,    24,    29,    49,
      47,     9,    23,    24,    25,    26,    27,    28,    44,    52,
      49,    16,    47,     8,     8,    24,    22,    24,    49,    47,
      22,    24,    44,     8,    24,    22,    24,    31,    33,    49,
      13,    48,    49,    22,    24,    48,    49,     8,    49,    48,
      49,    24,    48,    49,    24,    24,    49,    24,    49,    49,
      49
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
/* Line 1792 of yacc.c  */
#line 103 "chn_vsi.y"
    {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); YYACCEPT; 
		}
    break;

  case 3:
/* Line 1792 of yacc.c  */
#line 106 "chn_vsi.y"
    {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); free((yyvsp[(1) - (4)].pcval)); YYACCEPT;
		}
    break;

  case 4:
/* Line 1792 of yacc.c  */
#line 109 "chn_vsi.y"
    {
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); YYACCEPT; 
		}
    break;

  case 5:
/* Line 1792 of yacc.c  */
#line 112 "chn_vsi.y"
    { yyerror((yyvsp[(3) - (4)].Sval)); YYERROR; }
    break;

  case 6:
/* Line 1792 of yacc.c  */
#line 113 "chn_vsi.y"
    { yyerror((yyvsp[(2) - (3)].Sval)); YYERROR; }
    break;

  case 7:
/* Line 1792 of yacc.c  */
#line 114 "chn_vsi.y"
    { yyerror((yyvsp[(1) - (2)].Sval)); YYERROR; }
    break;

  case 8:
/* Line 1792 of yacc.c  */
#line 115 "chn_vsi.y"
    { yyerror("parse error"); YYERROR; }
    break;

  case 9:
/* Line 1792 of yacc.c  */
#line 118 "chn_vsi.y"
    { strcpy((yyval.Sval),(yyvsp[(2) - (2)].Sval)); }
    break;

  case 10:
/* Line 1792 of yacc.c  */
#line 119 "chn_vsi.y"
    { strcpy((yyval.Sval),(yyvsp[(1) - (1)].Sval)); }
    break;

  case 11:
/* Line 1792 of yacc.c  */
#line 122 "chn_vsi.y"
    { 
			(yyval.Lval)=(yyvsp[(2) - (3)].Lval); TAIL_NODE->next=(yyvsp[(2) - (3)].Lval); TAIL_NODE=(yyvsp[(2) - (3)].Lval); TAIL_NODE->comment=0;
		}
    break;

  case 12:
/* Line 1792 of yacc.c  */
#line 125 "chn_vsi.y"
    { 
			(yyval.Lval)=(yyvsp[(2) - (3)].Lval); TAIL_NODE->next=(yyvsp[(2) - (3)].Lval); TAIL_NODE=(yyvsp[(2) - (3)].Lval); TAIL_NODE->comment=0;
		}
    break;

  case 13:
/* Line 1792 of yacc.c  */
#line 128 "chn_vsi.y"
    { 
			(yyval.Lval)=(yyvsp[(2) - (4)].Lval); TAIL_NODE->next=(yyvsp[(2) - (4)].Lval); TAIL_NODE=(yyvsp[(2) - (4)].Lval); TAIL_NODE->comment=(yyvsp[(3) - (4)].pcval);
		}
    break;

  case 14:
/* Line 1792 of yacc.c  */
#line 131 "chn_vsi.y"
    { TAIL_NODE=HEAD_NODE=0; }
    break;

  case 15:
/* Line 1792 of yacc.c  */
#line 132 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(3) - (4)].Lval); TAIL_NODE=HEAD_NODE=(yyvsp[(3) - (4)].Lval); (yyval.Lval)->comment=0; }
    break;

  case 16:
/* Line 1792 of yacc.c  */
#line 133 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(3) - (4)].Lval); TAIL_NODE=HEAD_NODE=(yyvsp[(3) - (4)].Lval); (yyval.Lval)->comment=0; }
    break;

  case 17:
/* Line 1792 of yacc.c  */
#line 134 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(3) - (5)].Lval); TAIL_NODE=HEAD_NODE=(yyvsp[(3) - (5)].Lval); (yyval.Lval)->comment=(yyvsp[(4) - (5)].pcval); }
    break;

  case 18:
/* Line 1792 of yacc.c  */
#line 137 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 19:
/* Line 1792 of yacc.c  */
#line 138 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (2)].Lval); (yyvsp[(1) - (2)].Lval)->item.residue.width=(yyvsp[(2) - (2)].ival); }
    break;

  case 20:
/* Line 1792 of yacc.c  */
#line 139 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 21:
/* Line 1792 of yacc.c  */
#line 140 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (2)].Lval); (yyvsp[(1) - (2)].Lval)->item.molecule.width=(yyvsp[(2) - (2)].ival); }
    break;

  case 22:
/* Line 1792 of yacc.c  */
#line 141 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 23:
/* Line 1792 of yacc.c  */
#line 142 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 24:
/* Line 1792 of yacc.c  */
#line 143 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (2)].Lval); (yyvsp[(1) - (2)].Lval)->item.trace.width=(yyvsp[(2) - (2)].ival); }
    break;

  case 25:
/* Line 1792 of yacc.c  */
#line 144 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 26:
/* Line 1792 of yacc.c  */
#line 147 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNode("clear.",'C'); }
    break;

  case 27:
/* Line 1792 of yacc.c  */
#line 148 "chn_vsi.y"
    { (yyval.Lval)=(yyvsp[(1) - (1)].Lval); }
    break;

  case 28:
/* Line 1792 of yacc.c  */
#line 150 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNode((yyvsp[(1) - (1)].pcval),'E'); }
    break;

  case 30:
/* Line 1792 of yacc.c  */
#line 154 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeView((yyvsp[(1) - (1)].Vval)); }
    break;

  case 31:
/* Line 1792 of yacc.c  */
#line 157 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeTrace((yyvsp[(1) - (5)].ival),(yyvsp[(3) - (5)].ival),0,(yyvsp[(4) - (5)].cval),(yyvsp[(5) - (5)].cval),TRACE_WIDTH); }
    break;

  case 32:
/* Line 1792 of yacc.c  */
#line 159 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeTrace((yyvsp[(1) - (6)].ival),(yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].cval),(yyvsp[(5) - (6)].cval),(yyvsp[(6) - (6)].cval),TRACE_WIDTH); }
    break;

  case 33:
/* Line 1792 of yacc.c  */
#line 162 "chn_vsi.y"
    { (yyval.cval)='^'; }
    break;

  case 34:
/* Line 1792 of yacc.c  */
#line 163 "chn_vsi.y"
    { (yyval.cval)='.'; }
    break;

  case 35:
/* Line 1792 of yacc.c  */
#line 164 "chn_vsi.y"
    { (yyval.cval)='+'; }
    break;

  case 36:
/* Line 1792 of yacc.c  */
#line 165 "chn_vsi.y"
    { (yyval.cval)='~'; }
    break;

  case 37:
/* Line 1792 of yacc.c  */
#line 166 "chn_vsi.y"
    { (yyval.cval)=':'; }
    break;

  case 38:
/* Line 1792 of yacc.c  */
#line 167 "chn_vsi.y"
    { (yyval.cval)='@'; }
    break;

  case 39:
/* Line 1792 of yacc.c  */
#line 170 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (4)].cval),(yyvsp[(2) - (4)].ival),KEY_CHAIN,0,0,(yyvsp[(4) - (4)].Cval)[0],(yyvsp[(4) - (4)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 40:
/* Line 1792 of yacc.c  */
#line 172 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (5)].cval),(yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].cval),0,0,(yyvsp[(5) - (5)].Cval)[0],(yyvsp[(5) - (5)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 41:
/* Line 1792 of yacc.c  */
#line 174 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (6)].cval),(yyvsp[(2) - (6)].ival),KEY_CHAIN,(yyvsp[(4) - (6)].pcval),0,(yyvsp[(6) - (6)].Cval)[0],(yyvsp[(6) - (6)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 42:
/* Line 1792 of yacc.c  */
#line 176 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (7)].cval),(yyvsp[(2) - (7)].ival),(yyvsp[(3) - (7)].cval),(yyvsp[(5) - (7)].pcval),0,(yyvsp[(7) - (7)].Cval)[0],(yyvsp[(7) - (7)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 43:
/* Line 1792 of yacc.c  */
#line 179 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (6)].cval),(yyvsp[(2) - (6)].ival),KEY_CHAIN,AllocString("c"),AllocString("o"),
						(yyvsp[(6) - (6)].Cval)[0],(yyvsp[(6) - (6)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 44:
/* Line 1792 of yacc.c  */
#line 182 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (7)].cval),(yyvsp[(2) - (7)].ival),(yyvsp[(3) - (7)].cval),AllocString("c"), AllocString("o"),
						(yyvsp[(7) - (7)].Cval)[0],(yyvsp[(7) - (7)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 45:
/* Line 1792 of yacc.c  */
#line 187 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (8)].cval),(yyvsp[(2) - (8)].ival),KEY_CHAIN,(yyvsp[(4) - (8)].pcval),(yyvsp[(6) - (8)].pcval),(yyvsp[(8) - (8)].Cval)[0],(yyvsp[(8) - (8)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 46:
/* Line 1792 of yacc.c  */
#line 190 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeRes((yyvsp[(1) - (9)].cval),(yyvsp[(2) - (9)].ival),(yyvsp[(3) - (9)].cval),(yyvsp[(5) - (9)].pcval),(yyvsp[(7) - (9)].pcval),(yyvsp[(9) - (9)].Cval)[0],(yyvsp[(9) - (9)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 47:
/* Line 1792 of yacc.c  */
#line 194 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (4)].pcval),(yyvsp[(2) - (4)].ival),0,0,0,(yyvsp[(4) - (4)].Cval)[0],(yyvsp[(4) - (4)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 48:
/* Line 1792 of yacc.c  */
#line 196 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (5)].pcval),(yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].cval),0,0,(yyvsp[(5) - (5)].Cval)[0],(yyvsp[(5) - (5)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 49:
/* Line 1792 of yacc.c  */
#line 198 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (6)].pcval),(yyvsp[(2) - (6)].ival),0,(yyvsp[(4) - (6)].pcval),0,(yyvsp[(6) - (6)].Cval)[0],(yyvsp[(6) - (6)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 50:
/* Line 1792 of yacc.c  */
#line 200 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (7)].pcval),(yyvsp[(2) - (7)].ival),(yyvsp[(3) - (7)].cval),(yyvsp[(5) - (7)].pcval),0,(yyvsp[(7) - (7)].Cval)[0],(yyvsp[(7) - (7)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 51:
/* Line 1792 of yacc.c  */
#line 202 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (8)].pcval),(yyvsp[(2) - (8)].ival),0,(yyvsp[(4) - (8)].pcval),(yyvsp[(6) - (8)].pcval),(yyvsp[(8) - (8)].Cval)[0],(yyvsp[(8) - (8)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 52:
/* Line 1792 of yacc.c  */
#line 204 "chn_vsi.y"
    { (yyval.Lval)=MakeItemNodeMol((yyvsp[(1) - (9)].pcval),(yyvsp[(2) - (9)].ival),(yyvsp[(3) - (9)].cval),(yyvsp[(5) - (9)].pcval),(yyvsp[(7) - (9)].pcval),(yyvsp[(9) - (9)].Cval)[0],(yyvsp[(9) - (9)].Cval)[1],THIN_WIRE_WIDTH); }
    break;

  case 53:
/* Line 1792 of yacc.c  */
#line 208 "chn_vsi.y"
    { (yyval.pcval)=AllocString((yyvsp[(1) - (1)].aval)); }
    break;

  case 54:
/* Line 1792 of yacc.c  */
#line 210 "chn_vsi.y"
    { (yyval.pcval)=AllocString((yyvsp[(1) - (1)].aval)); }
    break;

  case 55:
/* Line 1792 of yacc.c  */
#line 212 "chn_vsi.y"
    { (yyval.Cval)[0]=(yyvsp[(2) - (3)].cval); (yyval.Cval)[1]='T'; }
    break;

  case 56:
/* Line 1792 of yacc.c  */
#line 213 "chn_vsi.y"
    { (yyval.Cval)[0]=(yyvsp[(2) - (3)].cval); (yyval.Cval)[1]='S'; }
    break;

  case 57:
/* Line 1792 of yacc.c  */
#line 214 "chn_vsi.y"
    { (yyval.Cval)[0]=(yyvsp[(1) - (1)].cval); (yyval.Cval)[1]='F'; }
    break;

  case 58:
/* Line 1792 of yacc.c  */
#line 217 "chn_vsi.y"
    { char s[3]; s[0]=(yyvsp[(2) - (2)].cval); s[1]=0; strcpy((yyval.Sval),(yyvsp[(1) - (2)].Sval)); strcat((yyval.Sval),s); 
			if( strlen((yyval.Sval)) > 90) yyerror((yyval.Sval)); }
    break;

  case 59:
/* Line 1792 of yacc.c  */
#line 219 "chn_vsi.y"
    { (yyval.Sval)[0]=(yyvsp[(1) - (1)].cval); (yyval.Sval)[1]=0; (yyval.Sval)[2]=0; }
    break;

  case 60:
/* Line 1792 of yacc.c  */
#line 221 "chn_vsi.y"
    { PDB_File_Comment=0; }
    break;

  case 61:
/* Line 1792 of yacc.c  */
#line 222 "chn_vsi.y"
    { PDB_File_Comment=(yyvsp[(2) - (2)].pcval); }
    break;

  case 62:
/* Line 1792 of yacc.c  */
#line 225 "chn_vsi.y"
    { free((yyvsp[(2) - (3)].pcval)); }
    break;

  case 63:
/* Line 1792 of yacc.c  */
#line 226 "chn_vsi.y"
    { }
    break;

  case 64:
/* Line 1792 of yacc.c  */
#line 227 "chn_vsi.y"
    { }
    break;


/* Line 1792 of yacc.c  */
#line 1896 "chn_vsi.tab.c"
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2055 of yacc.c  */
#line 230 "chn_vsi.y"


int	yyerror(char *s)
{ 
	char	str[102];
	if(FileName == 0){
		fprintf(stderr,"File number = %d\n",KEY_FILE);
		print_error("Invalid file number!");
	} else {
		fprintf(stderr,"\nfatal error (line %d): '%s'.\n",PARSELINE,s);
		fprintf(stderr,"%s --> \"%s\"\n",s,fgets(str,100,yyin));
		print_error("Input rejected.\n");
	}
	return 1;
}

void    Close(FILE *fptr) {if(fptr != stderr && fptr != stdout) fclose(fptr); }

#define USAGE_START "Usage: chn_vsi infile.vsi <int> [options]\n\
Usage: chn_vsi infile.vsi <int> <int> [options] (for multiple input files).\n\
       -A<int_array>  sort and print the files designated in <int_array>\n\
       -a<int>     print a copy of the current file adding <int> to residue sites\n\
       -append     append short rasmol add-on to standard output file\n\
       -Append     append full rasmol add-on to standard output file\n\
       -f<int>     print a copy of the current file adding file <int> to items\n\
       -c          show cis vs trans proline\n\
       -d<real>    dmax in Angstroms for classical H-bonds (default: 2.7 Angstoms)\n\
       -C          sort by chain and print a copy of the current file\n\
       -D          look for pdb input file within current directory\n\
       -B          Burley method for aromatic-aromatic interactions\n\
       -pml        generate a PyMol script instead\n\
       -pml=<file> generate a PyMol script instead\n\
       -s          sort and print a copy of the current file\n\
       -S          sort, purge, and print a copy of the current file\n\
       -skip=<char> if skip='W' the skip weak H-bonds; if skip == 'A' skip all H-bonds\n\
       -T          Test mode (fill in H-bonds)\n\
       -out=<name> Specify name of the output file\n\
       -t          Print protein regions w/o backbone traces (default: skip such regions)\n\
       -v          verbose mode\n\
       -w          print a copy of the current file\n\
       -x          dummy\n\n"


// int	main(int argc,char *argv[])
int	ChainVSI(int argc,char *argv[],FILE *ifp, FILE *ofp)
{ 
	Int4	i,n,arg,add=0,file=0,*value;
	char	c,mode=' ',filename[200]; filename[0]=0;
	vsi_typ	node,*headnode=0,*tailnode=0,tmp_node;
	item_type itm;
	char	**FileNames=0,*KeyChains=0,**FileComment=0;
	FILE	*fp;
	BooLean	only_if_trace=TRUE;
	Int4	append_output=0;
	char	modeHB=0,*pml_file=0,str[200];
        float   HA_dmax = 2.7,dmax=3.6;
	int	rtn_number=0;
	Int4	SARP_ID=0;	// the identifier for key file within a SARP multiple VSI input files.
	Int4	arg_start=3;

	if(argc < 3) print_error(USAGE_START);
	if(sscanf(argv[2],"%d",&i)!=1) print_error(USAGE_START);
	KEY_FILE=i;
	if(i <1) print_error(USAGE_START);
#if 0
        for(arg = 0; arg < argc; arg++) fprintf(stderr,"%s ",argv[arg]); 
	fprintf(stderr,"\n");
#endif
	if(argc >= 4 && argv[3][0] != '-'){
		if(sscanf(argv[3],"%d",&SARP_ID)!=1) print_error(USAGE_START); arg_start=4; 
		// fprintf(stderr,"sarp id = %d\n",SARP_ID);
	}
        for(arg = arg_start; arg < argc; arg++){
          if(argv[arg][0] == '-'){
           switch(argv[arg][1]) {
             case 'c': CIS_TRANS_PRO=TRUE; break;
             case 'D': USE_CURRENT_DIR=TRUE; break;
             case 'a': 
		if(strcmp("-append",argv[arg]) == 0){
			append_output=1; 
			// print_error(" -append option not yet implemented (problems)");
		} else {
		  mode='w'; 
		  add=IntOption(argv[arg],'a',-5000,5000,USAGE_START); 
		} break;
             case 'f': mode='w'; file=IntOption(argv[arg],'f',1,5000,USAGE_START); break;
             case 'w': mode='w'; break;
             case 'd': HA_dmax=RealOption(argv[arg],'d',1.0,50,USAGE_START); break;
             case 's': 
		if(sscanf(argv[arg],"-skip=%c",&modeHB) != 1){
			if(argv[arg][2]==0) mode='s'; else print_error(USAGE_START);
		} break;
             case 'S': mode='S'; break;
             case 'C': mode='C'; break;
             case 'B': mode='B'; break;
             case 'T': mode='T'; break;
             case 't': only_if_trace=FALSE; break;
             case 'v': VERBOSE=TRUE; break;
             case 'p':  mode='p';
		if(sscanf(argv[arg],"-pml=%s",str) == 1){ 
			pml_file=AllocString(str);
		} else if(strcmp("-pml",argv[arg]) == 0){ 
			pml_file=0;
		} else print_error(USAGE_START);
		break;
             case 'o': mode=' ';
		if(sscanf(argv[arg],"-out=%s",&filename) != 1) print_error(USAGE_START);
		break;
             case 'A': 
	      if(strcmp("-Append",argv[arg]) == 0){
			append_output=2; 
	      } else {
		mode='A'; 
		Int4 tmp_value[1000];
		n=ParseIntegers(argv[arg]+2,tmp_value,USAGE_START);
		NEW(value,n+3, Int4);
		for(i=0; i<= n; i++) value[i]=tmp_value[i];
		NEW(headnode,n+3, vsi_typ);
		NEW(tailnode,n+3, vsi_typ);
		NEWP(FileNames,n+3, char);
		NEWP(FileComment,n+3, char);
		NEW(KeyChains,n+3, char);
	      } break;
             default : print_error(USAGE_START);
           }
	  }else print_error(USAGE_START);
	}
	if(mode != 'A'){
	   NUMBER_OF_PDB_FILES=0;
	   if(debug) fprintf(stderr,"=========== Parsing file %s ============\n",argv[1]);
	   if(ifp) yyin=ifp; else yyin = open_file(argv[1],"","r"); yyout=stderr; 
	   if(SARP_ID >0){	// then input contains multiple vsi files; go to starting file..
		char str[200];
		char	c; 
		Int4	i;
		do {
		  while((c=fgetc(yyin)) != '~'){
			if(c == EOF){
				fprintf(stderr,"Input does not contain a vsi file with ID %d\n",SARP_ID);
				print_error("Fatal error");
			}
		  } fgets(str,198,yyin);
		  if(sscanf(str,"$=%d.\n",&i)!=1) print_error(USAGE_START);
		} while(i != SARP_ID);
	   }
	   yyparse(); 
	   if(ifp==0) fclose(yyin); 
	   rtn_number=NUMBER_OF_PDB_FILES;
	   if(debug) fprintf(stderr,"=========== %d files parsed. ============\n",NUMBER_OF_PDB_FILES);
	}
	switch (mode){
	  case 'A':
	     file=0;
	     for(i=1; i <=n; i++){
		KEY_FILE=value[i]; 
	   	NUMBER_OF_PDB_FILES=0;
		yyin = open_file(argv[1],"","r"); yyparse(); fclose(yyin);
	   	rtn_number=NUMBER_OF_PDB_FILES;
		// fprintf(stderr,"value[%d]=%d (%c)\n",i,value[i],KEY_CHAIN);
		headnode[i]=HEAD_NODE;
		tailnode[i]=TAIL_NODE;
		FileNames[i]=FileName;
		FileComment[i]=PDB_File_Comment;
		KeyChains[i]=KEY_CHAIN;
		for(node=HEAD_NODE; node; node=node->next){
		    node->file=value[i];
		    itm=node->item;
		    if(node->type == 'R'){
			if(itm.residue.chain == 0) node->item.residue.chain=KEY_CHAIN; 
			// else fprintf(stderr,"res chain: %c\n",itm.residue.chain);
		    } else if(node->type == 'T'){
			if(itm.trace.chain == 0) node->item.trace.chain=KEY_CHAIN;
			// else fprintf(stderr,"trace chain: %c\n",itm.residue.chain);
		    }
		}
		HEAD_NODE=TAIL_NODE=0; FileName=0; PARSELINE=1;
	     }
	     HEAD_NODE=headnode[1];
	     for(i=1; i < n; i++){ tailnode[i]->next=headnode[i+1]; }
	     KEY_FILE=0;
	     HEAD_NODE=SortItems(HEAD_NODE);
	     mode='C';
	  case 'C':
	  case 'S': case 's':
	     if(mode=='S') HEAD_NODE=SortItems(HEAD_NODE,TRUE);
	     else if(mode=='C') HEAD_NODE=SortItemsByChains(HEAD_NODE);
	     else HEAD_NODE=SortItems(HEAD_NODE);
	  case 'w':
	     fp=open_file(argv[1],".crs","w");
	     if(FileNames){
	     	PrintItems(fp,HEAD_NODE,n,FileNames,FileComment,KeyChains,value,add,0,FALSE);
	     } else {
	     	PrintItems(fp,HEAD_NODE,FileName,PDB_File_Comment,KEY_CHAIN,KEY_FILE,
				add,file,FALSE);
	     }
	     Close(fp);
	     break;
	  case 'B':
	     {
//***************** Begin Burley Routine ****************************************
		char str[200],*str2,*DATADIR=0;
		BooLean *skip=0;

	    if(!USE_CURRENT_DIR){
        	DATADIR=getenv("CHN_RAS_DIR");
        	if(debug) fprintf(stderr,"CHN_RAS_DIR = %s\n",DATADIR);
		// if(DATADIR == 0) sprintf(str,"%s",FileName); else 
		sprintf(str,"%s/%s",DATADIR,FileName);
	    } else { sprintf(str,"%s",FileName); }
		str2=strstr(str,".pdb");
		*str2 = 0;
        	// fprintf(stderr,"pdb path = %s\n",str);
		pdb_typ P=MkPDB(str);
	        ErrorToDevNullPDB(P);

        Int4    C=GetChainNumberPDB(P,KEY_CHAIN);
        Int4    resStart=MinResPDB(C,P);
        Int4    resEnd=MaxResPDB(C,P),r;

		NEW(skip,resEnd+5,BooLean);
		for(r=0; r <= resEnd; r++) skip[r]=TRUE;
		for(node=HEAD_NODE; node; node=node->next){
		    node->file=KEY_FILE;
		    itm=node->item;
		    if(node->type == 'R'){
			if(itm.residue.chain == 0) node->item.residue.chain=KEY_CHAIN; 
			if(itm.residue.chain != KEY_CHAIN) continue;
			r=itm.residue.site;
			if(!(r >= resStart && r <= resEnd)){
			  fprintf(stderr,
			  "WARNING: Inconsistent input r=%d; resStart=%d; resEnd=%d\n",
					r,resStart,resEnd);
			  fprintf(stderr," (ignored)\n");
			} else skip[r]=FALSE;
			// else fprintf(stderr,"res chain: %c\n",itm.residue.chain);
		    }
		}
   {
        Int4    file_id=0;
        char    color=0;

	file_id=KEY_FILE;
        if(KEY_CHAIN==0) KEY_CHAIN='X';
#if 0
        PutAroAroInteractions(stdout,file_id,resStart,resEnd,HA_dmax,dmax,C,
                C,KEY_CHAIN,skip,color,P);
#endif
	free(skip);
   }
		NilPDB(P);
  }
//***************** End Burley Routine ****************************************
	     break;
//***************** Begin AutoFormating ***************************************
	  case 'T':
	{
	    char    str[200],*str2,*DATADIR=0;
	    Int4    N;
	    FILE    *fp;

	    if(!USE_CURRENT_DIR){
	      // if(FileName[0] == '/') DATADIR = AllocString(""); else DATADIR=getenv("CHN_RAS_DIR");
	      DATADIR=getenv("CHN_RAS_DIR");
	      if(debug) fprintf(stderr,"CHN_RAS_DIR = %s\n",DATADIR);
	      sprintf(str,"%s/%s",DATADIR,FileName);
	    } else sprintf(str,"%s",FileName);
#if 0
	    str2=strstr(str,".pdb"); *str2 = 0; 
	    // fprintf(stderr,"pdb path = %s\n",str);
	    pdb_typ P=MkPDB(str);
#else
	    pdb_typ P=MakePDB(str);
#endif
	    ErrorToDevNullPDB(P);
		
	    if(ofp) fp=tmpfile(); else
	    if(!USE_CURRENT_DIR){ fp=open_file(argv[1],"_tmp.vsi","w"); }
	    else { sprintf(str,"_tmp%d.vsi",KEY_FILE); fp=open_file(argv[1],str,"w"); }

	    if(KEY_CHAIN==0) fprintf(fp,"File%d=%s:X   \n\n",KEY_FILE,FileName);
	    else fprintf(fp,"File%d=%s:%c   \n\n",KEY_FILE,FileName,KEY_CHAIN);
	    vsi_typ first_head,first_tail=0;
	    first_head=PutAutoFormatVSI(fp,KEY_FILE,HEAD_NODE,KEY_CHAIN,
						HA_dmax,dmax, only_if_trace,P,modeHB);
	    if(ofp==0) fclose(fp); else { rewind(fp); }

	    for(N=0,node=first_head; node; node=node->next){ first_tail=node; N++; }
	    if(N <= 1) fprintf(stderr,"WARNING: ");
	    // fprintf(stderr,"%d items parsed\n",N);

	    if(PDB_File_Comment) free(PDB_File_Comment); free(FileName); 
	    FileName=0; PDB_File=0; PDB_File_Comment=0; KEY_CHAIN=0; PARSELINE=1; 
	    HEAD_NODE=TAIL_NODE=0; 
	    if(ofp){ yyin=fp; } else if(!USE_CURRENT_DIR){ yyin = open_file(argv[1],"_tmp.vsi","r"); }
	    else {
	        if(debug) fprintf(stderr,"=========== Parsing file %s_tmp%d.vsi ============\n",argv[1],KEY_FILE);
		sprintf(str,"_tmp%d.vsi",KEY_FILE); yyin = open_file(argv[1],str,"r"); 
	    }
	    NUMBER_OF_PDB_FILES=0;
	    yyout=stderr; yyparse(); fclose(yyin);	// don't close input file here.. 
	    if(debug) fprintf(stderr,"=========== %d files parsed. ============\n",NUMBER_OF_PDB_FILES);

	    for(N=0,node=HEAD_NODE; node; node=node->next){ node->file=KEY_FILE; N++; }
	    first_tail->next = HEAD_NODE; 
	    // HEAD_NODE=SortItems(first_head);		// sort nodes only...
	    HEAD_NODE=SortItems(first_head,TRUE);  // both sort & purge the nodes
	    if(N <= 1) fprintf(stderr,"WARNING: ");

	    // fprintf(stderr,"%d items parsed\n",N);
	    if(ofp) fp=ofp;
	    else if(!USE_CURRENT_DIR){
		fp=open_file(argv[1],".crs","w"); 
	        if(debug) fprintf(stderr,"=========== Creating file %s.crs ...============\n",argv[1]);
	    // } else { sprintf(str,"%d.vsi",KEY_FILE); fp=open_file(argv[1],str,"w"); }
	    } else { sprintf(str,".crs"); fp=open_file(argv[1],str,"w"); }
	    PrintItems(fp,HEAD_NODE,FileName,PDB_File_Comment,KEY_CHAIN,KEY_FILE,add,file,FALSE);
	    if(ofp == 0) Close(fp);

	    if(USE_CURRENT_DIR){
	        if(debug) fprintf(stderr,"======= Interaction network file '%s%d.vsi' created. =======\n",
			argv[1],KEY_FILE);
	    }
	    for(node=HEAD_NODE; node; ){ tmp_node=node->next; FreeNode(node);  node=tmp_node; }
	    if(PDB_File_Comment) free(PDB_File_Comment); free(FileName); NilPDB(P);
	}
	     break;
//******************* END AutoFormating ***************************************
//*****************************************************************************
	  case 'r':
	   fp=open_file(filename,"","w");
	   PrintItemsRasmol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
#if 1
	   if(append_output){
          	fprintf(fp,"select hetero and within(4.0,*%c)\n",KEY_CHAIN);
          	fprintf(fp,"wireframe 50\n");
          	fprintf(fp,"select water\ncolor cpk\n");
          	fprintf(fp,"select (ions,Mg,Zn,Mn,Fe,Co) and within(4.0,*%c)\n",
                                                                KEY_CHAIN);
          	fprintf(fp,"color cpk\nspacefill 250\n");
		if(append_output == 2){
          	  fprintf(fp,"select (protein or water or (ions,Mg,Zn,Mn,Fe,Co)) ");
          	  fprintf(fp,"and not within(4.0,*%c)\n",KEY_CHAIN);
          	  fprintf(fp,"spacefill off\nwireframe off\n");
		}
	   }
#endif
	   Close(fp);
	   fprintf(stderr,"======= Rasmol script file %s created. =======\n",filename);
	   break;
	  case ' ':
	  case 'p':
	  default:
	   if(filename[0]==0) sprintf(filename,"%s",argv[1]);
	   if(mode=='p'){
	     if(pml_file){
		fp=open_file(pml_file,".pml","w"); 
		fprintf(stderr,"=========== Creating file %s.pml ============\n",pml_file);
		free(pml_file); pml_file=0;
	     } else {
		fp=stdout;
		fprintf(stderr,"=========== Creating file %s.pml ============\n",argv[1]);
	     }
	     PrintItemsPyMol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
	   } else {	// rasmol script instead...
	     fp=open_file(filename,".ras","w");
	     fprintf(stderr,"=========== Creating file %s.ras ============\n",argv[1]);
	     PrintItemsRasmol(fp,!USE_CURRENT_DIR,HEAD_NODE,FileName,KEY_CHAIN,
					CIS_TRANS_PRO,WIRE_WIDTH,SPACEFILL);
	   } Close(fp); 
	   if(PDB_File_Comment) free(PDB_File_Comment); free(FileName);
	   for(node=HEAD_NODE; node; ){ tmp_node=node->next; FreeNode(node);  node=tmp_node; }
	   break;
	}
	return rtn_number;
}

