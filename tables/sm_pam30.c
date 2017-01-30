/*  $Id: sm_pam30.c 90506 2006-09-25 19:30:59Z madden $
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

/** Entries for the PAM30 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Pam30PSM[25 * 25] = {
 /*       A,,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M , F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */  
/*A*/ 4,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*R*/ -30,5,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*N*/ -30,-30,6,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*D*/ -30,-30,-30,6,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*C*/ -30,-30,-30,-30,9,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*Q*/ -30,-30,-30,-30,-30,5,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*E*/ -30,-30,-30,-30,-30,-30,5,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*G*/ -30,-30,-30,-30,-30,-30,-30,6,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*H*/ -30,-30,-30,-30,-30,-30,-30,-30,8,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*I*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*L*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*K*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,5,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*M*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,5,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*F*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,6,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*P*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,7,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*S*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30,-30,-30,-30,-30,-30,-30,-30, 
/*T*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,5,-30,-30,-30,-30,-30,-30,-30,-30, 
/*W*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,11,-30,-30,-30,-30,-30,-30,-30, 
/*Y*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,7,-30,-30,-30,-30,-30,-30, 
/*V*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30,-30,-30,-30, 
/*B*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30,-30,-30, 
/*J*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,3,-30,-30,-30, 
/*Z*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,4,-30,-30, 
/*X*/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-1,-30, 
/***/ -30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,1
};
const SNCBIPackedScoreMatrix NCBISM_Pam30 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Pam30PSM,
    -17
};
