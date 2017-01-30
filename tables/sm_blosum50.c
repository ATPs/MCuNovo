/*  $Id: sm_blosum50.c 90507 2006-09-25 19:31:51Z madden $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Author:  Aaron Ucko, Mike Gertz
*
* File Description:
*   Protein alignment score matrices; shared between the two toolkits.
*
* ===========================================================================
*/

#include <util/tables/raw_scoremat.h>

/** Entries for the BLOSUM50 matrix at a scale of ln(2)/3.0. */

static const TNCBIScore s_Blosum50PSM[25 * 25] = {
  /*       A  ,   R,   N,   D,   C,   Q,   E,   G,   H,   I,   L,   K,   M,   F,   P,   S,   T,   W,   Y,   V,   B,   J,   Z,   X,   *        */ 
/*A*/  4, -5, -10, -10, 0, -5, -5, 0, -10, -5, -5, -5, -5, -10, -5, 1, 0, -15, -10, 0, -10, -5, -5, -5, -20,
/*R*/  -5, 5, 0, -10, -15, 1, 0, -10, 0, -15, -10, 2, -5, -15, -10, -5, -5, -15, -10, -15, -5, -10, 0, -5, -20,
/*N*/  -10, 0, 6, 1, -15, 0, 0, 0, 1, -15, -15, 0, -10, -15, -10, 1, 0, -20, -10, -15, 4, -15, 0, -5, -20,
/*D*/  -10, -10, 1, 6, -15, 0, 2, -5, -5, -15, -20, -5, -15, -15, -5, 0, -5, -20, -15, -15, 4, -15, 1, -5, -20,
/*C*/  0, -15, -15, -15, 9, -15, -20, -15, -15, -5, -5, -15, -5, -10, -15, -5, -5, -10, -10, -5, -15, -5, -15, -5, -20,
/*Q*/  -5, 1, 0, 0, -15, 5, 2, -10, 0, -15, -10, 1, 0, -15, -5, 0, -5, -10, -5, -10, 0, -10, 4, -5, -20,
/*E*/  -5, 0, 0, 2, -20, 2, 5, -10, 0, -15, -15, 1, -10, -15, -5, 0, -5, -15, -10, -10, 1, -15, 4, -5, -20,
/*G*/  0, -10, 0, -5, -15, -10, -10, 6, -10, -20, -20, -10, -15, -15, -10, 0, -10, -10, -15, -15, -5, -20, -10, -5, -20,
/*H*/  -10, 0, 1, -5, -15, 0, 0, -10, 8, -15, -15, -5, -10, -5, -10, -5, -10, -10, 2, -15, 0, -15, 0, -5, -20,
/*I*/  -5, -15, -15, -15, -5, -15, -15, -20, -15, 4, 2, -15, 1, 0, -15, -10, -5, -15, -5, 3, -15, 3, -15, -5, -20,
/*L*/  -5, -10, -15, -20, -5, -10, -15, -20, -15, 2, 4, -10, 2, 0, -15, -10, -5, -10, -5, 1, -20, 3, -15, -5, -20,
/*K*/  -5, 2, 0, -5, -15, 1, 1, -10, -5, -15, -10, 5, -5, -15, -5, 0, -5, -15, -10, -10, 0, -15, 1, -5, -20,
/*M*/  -5, -5, -10, -15, -5, 0, -10, -15, -10, 1, 2, -5, 5, 0, -10, -5, -5, -5, -5, 1, -15, 2, -5, -5, -20,
/*F*/  -10, -15, -15, -15, -10, -15, -15, -15, -5, 0, 0, -15, 0, 6, -20, -10, -10, 1, 3, -5, -15, 0, -15, -5, -20,
/*P*/  -5, -10, -10, -5, -15, -5, -5, -10, -10, -15, -15, -5, -10, -20, 7, -5, -5, -20, -15, -10, -10, -15, -5, -5, -20,
/*S*/  1, -5, 1, 0, -5, 0, 0, 0, -5, -10, -10, 0, -5, -10, -5, 4, 1, -15, -10, -10, 0, -10, 0, -5, -20,
/*T*/  0, -5, 0, -5, -5, -5, -5, -10, -10, -5, -5, -5, -5, -10, -5, 1, 5, -10, -10, 0, -5, -5, -5, -5, -20,
/*W*/  -15, -15, -20, -20, -10, -10, -15, -10, -10, -15, -10, -15, -5, 1, -20, -15, -10, 11, 2, -15, -20, -10, -10, -5, -20,
/*Y*/  -10, -10, -10, -15, -10, -5, -10, -15, 2, -5, -5, -10, -5, 3, -15, -10, -10, 2, 7, -5, -15, -5, -10, -5, -20,
/*V*/  0, -15, -15, -15, -5, -10, -10, -15, -15, 3, 1, -10, 1, -5, -10, -10, 0, -15, -5, 4, -15, 2, -10, -5, -20,
/*B*/  -10, -5, 4, 4, -15, 0, 1, -5, 0, -15, -20, 0, -15, -15, -10, 0, -5, -20, -15, -15, 4, -15, 0, -5, -20,
/*J*/  -5, -10, -15, -15, -5, -10, -15, -20, -15, 3, 3, -15, 2, 0, -15, -10, -5, -10, -5, 2, -15, 3, -15, -5, -20,
/*Z*/  -5, 0, 0, 1, -15, 4, 4, -10, 0, -15, -15, 1, -5, -15, -5, 0, -5, -10, -10, -10, 0, -15, 4, -5, -20,
/*X*/  -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -1, -20,
/***/  -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, 1
};
const SNCBIPackedScoreMatrix NCBISM_Blosum50 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum50PSM,
    -5
};

