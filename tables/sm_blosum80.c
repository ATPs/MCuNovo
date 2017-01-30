/*  $Id: sm_blosum80.c 90506 2006-09-25 19:30:59Z madden $
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

/** Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Blosum80PSM[25 * 25] = {
/*       A, ,   R,   N,   D,   C,   Q,   E,   G,   H,   I,   L,   K,   M ,  F,   P,   S,   T,   W,   Y,   V,   B,   J,   Z,   X,   *        */
/*A*/  4, -3, -6, -6, 0, -3, -3, 0, -6, -3, -3, -3, -3, -6, -3, 1, 0, -9, -6, 0, -6, -3, -3, -3, -12,
/*R*/  -3, 5, 0, -6, -9, 1, 0, -6, 0, -9, -6, 2, -3, -9, -6, -3, -3, -9, -6, -9, -3, -6, 0, -3, -12,
/*N*/  -6, 0, 6, 1, -9, 0, 0, 0, 1, -9, -9, 0, -6, -9, -6, 1, 0, -12, -6, -9, 4, -9, 0, -3, -12,
/*D*/  -6, -6, 1, 6, -9, 0, 2, -3, -3, -9, -12, -3, -9, -9, -3, 0, -3, -12, -9, -9, 4, -9, 1, -3, -12,
/*C*/  0, -9, -9, -9, 9, -9, -12, -9, -9, -3, -3, -9, -3, -6, -9, -3, -3, -6, -6, -3, -9, -3, -9, -3, -12,
/*Q*/  -3, 1, 0, 0, -9, 5, 2, -6, 0, -9, -6, 1, 0, -9, -3, 0, -3, -6, -3, -6, 0, -6, 4, -3, -12,
/*E*/  -3, 0, 0, 2, -12, 2, 5, -6, 0, -9, -9, 1, -6, -9, -3, 0, -3, -9, -6, -6, 1, -9, 4, -3, -12,
/*G*/  0, -6, 0, -3, -9, -6, -6, 6, -6, -12, -12, -6, -9, -9, -6, 0, -6, -6, -9, -9, -3, -12, -6, -3, -12,
/*H*/  -6, 0, 1, -3, -9, 0, 0, -6, 8, -9, -9, -3, -6, -3, -6, -3, -6, -6, 2, -9, 0, -9, 0, -3, -12,
/*I*/  -3, -9, -9, -9, -3, -9, -9, -12, -9, 4, 2, -9, 1, 0, -9, -6, -3, -9, -3, 3, -9, 3, -9, -3, -12,
/*L*/  -3, -6, -9, -12, -3, -6, -9, -12, -9, 2, 4, -6, 2, 0, -9, -6, -3, -6, -3, 1, -12, 3, -9, -3, -12,
/*K*/  -3, 2, 0, -3, -9, 1, 1, -6, -3, -9, -6, 5, -3, -9, -3, 0, -3, -9, -6, -6, 0, -9, 1, -3, -12,
/*M*/  -3, -3, -6, -9, -3, 0, -6, -9, -6, 1, 2, -3, 5, 0, -6, -3, -3, -3, -3, 1, -9, 2, -3, -3, -12,
/*F*/  -6, -9, -9, -9, -6, -9, -9, -9, -3, 0, 0, -9, 0, 6, -12, -6, -6, 1, 3, -3, -9, 0, -9, -3, -12,
/*P*/  -3, -6, -6, -3, -9, -3, -3, -6, -6, -9, -9, -3, -6, -12, 7, -3, -3, -12, -9, -6, -6, -9, -3, -3, -12,
/*S*/  1, -3, 1, 0, -3, 0, 0, 0, -3, -6, -6, 0, -3, -6, -3, 4, 1, -9, -6, -6, 0, -6, 0, -3, -12,
/*T*/  0, -3, 0, -3, -3, -3, -3, -6, -6, -3, -3, -3, -3, -6, -3, 1, 5, -6, -6, 0, -3, -3, -3, -3, -12,
/*W*/  -9, -9, -12, -12, -6, -6, -9, -6, -6, -9, -6, -9, -3, 1, -12, -9, -6, 11, 2, -9, -12, -6, -6, -3, -12,
/*Y*/  -6, -6, -6, -9, -6, -3, -6, -9, 2, -3, -3, -6, -3, 3, -9, -6, -6, 2, 7, -3, -9, -3, -6, -3, -12,
/*V*/  0, -9, -9, -9, -3, -6, -6, -9, -9, 3, 1, -6, 1, -3, -6, -6, 0, -9, -3, 4, -9, 2, -6, -3, -12,
/*B*/  -6, -3, 4, 4, -9, 0, 1, -3, 0, -9, -12, 0, -9, -9, -6, 0, -3, -12, -9, -9, 4, -9, 0, -3, -12,
/*J*/  -3, -6, -9, -9, -3, -6, -9, -12, -9, 3, 3, -9, 2, 0, -9, -6, -3, -6, -3, 2, -9, 3, -9, -3, -12,
/*Z*/  -3, 0, 0, 1, -9, 4, 4, -6, 0, -9, -9, 1, -3, -9, -3, 0, -3, -6, -6, -6, 0, -9, 4, -3, -12,
/*X*/  -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -1, -12,
/***/  -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, 1
};
const SNCBIPackedScoreMatrix NCBISM_Blosum80 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum80PSM,
    -6
};

