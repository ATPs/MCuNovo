/*  $Id: sm_blosum90.c 90507 2006-09-25 19:31:51Z madden $
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

/** Entries for the BLOSUM90 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Blosum90PSM[25 * 25] = {
/*       A, ,   R,   N,   D,   C,   Q,   E,   G,   H,   I,   L,   K,   M ,  F,   P,   S,   T,   W,   Y,   V,   B,   J,   Z,   X,   *        */ 
/*A*/  4, -6, -7, -7, -5, -6, -6, -5, -7, -6, -6, -6, -6, -7, -6, -4, -5, -8, -7, -5, -7, -6, -6, -6, -9,
/*R*/  -6, 5, -5, -7, -8, -4, -5, -7, -5, -8, -7, -3, -6, -8, -7, -6, -6, -8, -7, -8, -6, -7, -5, -6, -9,
/*N*/  -7, -5, 6, -4, -8, -5, -5, -5, -4, -8, -8, -5, -7, -8, -7, -4, -5, -9, -7, -8, -1, -8, -5, -6, -9,
/*D*/  -7, -7, -4, 6, -8, -5, -3, -6, -6, -8, -9, -6, -8, -8, -6, -5, -6, -9, -8, -8, -1, -8, -4, -6, -9,
/*C*/  -5, -8, -8, -8, 9, -8, -9, -8, -8, -6, -6, -8, -6, -7, -8, -6, -6, -7, -7, -6, -8, -6, -8, -6, -9,
/*Q*/  -6, -4, -5, -5, -8, 5, -3, -7, -5, -8, -7, -4, -5, -8, -6, -5, -6, -7, -6, -7, -5, -7, -1, -6, -9,
/*E*/  -6, -5, -5, -3, -9, -3, 5, -7, -5, -8, -8, -4, -7, -8, -6, -5, -6, -8, -7, -7, -4, -8, -1, -6, -9,
/*G*/  -5, -7, -5, -6, -8, -7, -7, 6, -7, -9, -9, -7, -8, -8, -7, -5, -7, -7, -8, -8, -6, -9, -7, -6, -9,
/*H*/  -7, -5, -4, -6, -8, -5, -5, -7, 8, -8, -8, -6, -7, -6, -7, -6, -7, -7, -3, -8, -5, -8, -5, -6, -9,
/*I*/  -6, -8, -8, -8, -6, -8, -8, -9, -8, 4, -3, -8, -4, -5, -8, -7, -6, -8, -6, -2, -8, -2, -8, -6, -9,
/*L*/  -6, -7, -8, -9, -6, -7, -8, -9, -8, -3, 4, -7, -3, -5, -8, -7, -6, -7, -6, -4, -9, -2, -8, -6, -9,
/*K*/  -6, -3, -5, -6, -8, -4, -4, -7, -6, -8, -7, 5, -6, -8, -6, -5, -6, -8, -7, -7, -5, -8, -4, -6, -9,
/*M*/  -6, -6, -7, -8, -6, -5, -7, -8, -7, -4, -3, -6, 5, -5, -7, -6, -6, -6, -6, -4, -8, -3, -6, -6, -9,
/*F*/  -7, -8, -8, -8, -7, -8, -8, -8, -6, -5, -5, -8, -5, 6, -9, -7, -7, -4, -2, -6, -8, -5, -8, -6, -9,
/*P*/  -6, -7, -7, -6, -8, -6, -6, -7, -7, -8, -8, -6, -7, -9, 7, -6, -6, -9, -8, -7, -7, -8, -6, -6, -9,
/*S*/  -4, -6, -4, -5, -6, -5, -5, -5, -6, -7, -7, -5, -6, -7, -6, 4, -4, -8, -7, -7, -5, -7, -5, -6, -9,
/*T*/  -5, -6, -5, -6, -6, -6, -6, -7, -7, -6, -6, -6, -6, -7, -6, -4, 5, -7, -7, -5, -6, -6, -6, -6, -9,
/*W*/  -8, -8, -9, -9, -7, -7, -8, -7, -7, -8, -7, -8, -6, -4, -9, -8, -7, 11, -3, -8, -9, -7, -7, -6, -9,
/*Y*/  -7, -7, -7, -8, -7, -6, -7, -8, -3, -6, -6, -7, -6, -2, -8, -7, -7, -3, 7, -6, -8, -6, -7, -6, -9,
/*V*/  -5, -8, -8, -8, -6, -7, -7, -8, -8, -2, -4, -7, -4, -6, -7, -7, -5, -8, -6, 4, -8, -3, -7, -6, -9,
/*B*/  -7, -6, -1, -1, -8, -5, -4, -6, -5, -8, -9, -5, -8, -8, -7, -5, -6, -9, -8, -8, 4, -8, -5, -6, -9,
/*J*/  -6, -7, -8, -8, -6, -7, -8, -9, -8, -2, -2, -8, -3, -5, -8, -7, -6, -7, -6, -3, -8, 3, -8, -6, -9,
/*Z*/  -6, -5, -5, -4, -8, -1, -1, -7, -5, -8, -8, -4, -6, -8, -6, -5, -6, -7, -7, -7, -5, -8, 4, -6, -9,
/*X*/  -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -1, -9,
/***/  -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, 1
};
const SNCBIPackedScoreMatrix NCBISM_Blosum90 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum90PSM,
    -6
};

