#include "lulesh.h"

// If no MPI, then this whole file is stubbed out
#if USE_MPI

#include <mpi.h>
#include <string.h>
#include <iostream>

/* Comm Routines */

int cache_coherence_pad_real = CACHE_COHERENCE_PAD_REAL;

void RecvCommon(Real_t *data, int recvCount, MPI_Datatype &baseType, 
        int fromRank, int msgType, MPI_Comm comm, MPI_Request *req, int recvKind = 0)
{
    /* contiguous memory */
    //std::cout << "Waiting for irecv to complete." << std::endl;
    //std::cerr << "Receiving " << recvCount << " reals starting on " << data << " from " << fromRank << " with tag " << msgType << std::endl;
    MPI_Irecv(data, recvCount, baseType, fromRank, msgType, comm, req);
}

void SendCommon(Real_t *data, int sendCount, MPI_Datatype &baseType, 
        int destRank, int msgType, MPI_Comm comm, MPI_Request *req, int recvKind = 0)
{
    //std::cout << "Waiting for isend to complete." << std::endl;
    //std::cerr << "Sending " << sendCount << " reals starting on " << data << " to " << destRank << " with tag " << msgType << std::endl;
    //for (int i = 0; i < sendCount; i++) {
    //    std::cout << "  " << data[i];
    //}
    //std::cout << std::endl;
    MPI_Isend(data, sendCount, baseType, destRank, msgType, comm, req);
}

/* doRecv flag only works with regular block structure */
void CommRecv(Domain& domain, int msgType, Index_t xferFields,
        Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly, int it) 
{

    if (domain.numRanks() == 1)
        return ;

    /* post recieve buffers for all incoming messages */
    int myRank ;
    Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
    Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
    Index_t pmsg = 0 ; /* plane comm msg */
    Index_t emsg = 0 ; /* edge comm msg */
    Index_t cmsg = 0 ; /* corner comm msg */
    MPI_Datatype baseType = ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE) ;
    bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

    /* assume communication to 6 neighbors by default */
    rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;

    if (domain.rowLoc() == 0) {
        rowMin = false ;
    }
    if (domain.rowLoc() == (domain.tp()-1)) {
        rowMax = false ;
    }
    if (domain.colLoc() == 0) {
        colMin = false ;
    }
    if (domain.colLoc() == (domain.tp()-1)) {
        colMax = false ;
    }
    if (domain.planeLoc() == 0) {
        planeMin = false ;
    }
    if (domain.planeLoc() == (domain.tp()-1)) {
        planeMax = false ;
    }

    for (Index_t i=0; i<26; ++i) {
        domain.recvRequest[26*((msgType/8192)-1)+(i)] = MPI_REQUEST_NULL ;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

    /* post receives */

    /* receive data from neighboring domain faces */
    if (planeMin && doRecv) {
        /* contiguous memory */
        int fromRank = myRank - domain.tp()*domain.tp() ;
        int recvCount = dx * dy * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 1) ;
        ++pmsg ;
    }
    if (planeMax) {
        /* contiguous memory */
        int fromRank = myRank + domain.tp()*domain.tp() ;
        int recvCount = dx * dy * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 2) ;
        ++pmsg ;
    }
    if (rowMin && doRecv) {
        /* semi-contiguous memory */
        int fromRank = myRank - domain.tp() ;
        int recvCount = dx * dz * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 3) ;
        ++pmsg ;
    }
    if (rowMax) {
        /* semi-contiguous memory */
        int fromRank = myRank + domain.tp() ;
        int recvCount = dx * dz * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 4) ;
        ++pmsg ;
    }
    if (colMin && doRecv) {
        /* scattered memory */
        int fromRank = myRank - 1 ;
        int recvCount = dy * dz * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 5) ;
        ++pmsg ;
    }
    if (colMax) {
        /* scattered memory */
        int fromRank = myRank + 1 ;
        int recvCount = dy * dz * xferFields ;
        RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm],
                recvCount, baseType, fromRank, msgType,
                MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg)], 6) ;
        ++pmsg ;
    }

    if (!planeOnly) {
        /* receive data from domains connected only by an edge */
        if (rowMin && colMin && doRecv) {
            int fromRank = myRank - domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dz * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 7) ;
            ++emsg ;
        }

        if (rowMin && planeMin && doRecv) {
            int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dx * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 8) ;
            ++emsg ;
        }

        if (colMin && planeMin && doRecv) {
            int fromRank = myRank - domain.tp()*domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dy * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 9) ;
            ++emsg ;
        }

        if (rowMax && colMax) {
            int fromRank = myRank + domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dz * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 10) ;
            ++emsg ;
        }

        if (rowMax && planeMax) {
            int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dx * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 11) ;
            ++emsg ;
        }

        if (colMax && planeMax) {
            int fromRank = myRank + domain.tp()*domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dy * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 12) ;
            ++emsg ;
        }

        if (rowMax && colMin) {
            int fromRank = myRank + domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dz * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 13) ;
            ++emsg ;
        }

        if (rowMin && planeMax) {
            int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dx * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 14) ;
            ++emsg ;
        }

        if (colMin && planeMax) {
            int fromRank = myRank + domain.tp()*domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dy * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 15) ;
            ++emsg ;
        }

        if (rowMin && colMax && doRecv) {
            int fromRank = myRank - domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dz * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 16) ;
            ++emsg ;
        }

        if (rowMax && planeMin && doRecv) {
            int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dx * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 17) ;
            ++emsg ;
        }

        if (colMax && planeMin && doRecv) {
            int fromRank = myRank - domain.tp()*domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm],
                    dy * xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 18) ;
            ++emsg ;
        }

        /* receive data from domains connected only by a corner */
        if (rowMin && colMin && planeMin && doRecv) {
            /* corner at domain logical coord (0, 0, 0) */
            int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 19) ;
            ++cmsg ;
        }
        if (rowMin && colMin && planeMax) {
            /* corner at domain logical coord (0, 0, 1) */
            int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 20) ;
            ++cmsg ;
        }
        if (rowMin && colMax && planeMin && doRecv) {
            /* corner at domain logical coord (1, 0, 0) */
            int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 21) ;
            ++cmsg ;
        }
        if (rowMin && colMax && planeMax) {
            /* corner at domain logical coord (1, 0, 1) */
            int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 22) ;
            ++cmsg ;
        }
        if (rowMax && colMin && planeMin && doRecv) {
            /* corner at domain logical coord (0, 1, 0) */
            int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 23) ;
            ++cmsg ;
        }
        if (rowMax && colMin && planeMax) {
            /* corner at domain logical coord (0, 1, 1) */
            int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 24) ;
            ++cmsg ;
        }
        if (rowMax && colMax && planeMin && doRecv) {
            /* corner at domain logical coord (1, 1, 0) */
            int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 25) ;
            ++cmsg ;
        }
        if (rowMax && colMax && planeMax) {
            /* corner at domain logical coord (1, 1, 1) */
            int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;
            RecvCommon(&domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm +
                    emsg * maxEdgeComm +
                    cmsg * CACHE_COHERENCE_PAD_REAL],
                    xferFields, baseType, fromRank, msgType,
                    MPI_COMM_WORLD, &domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 26) ;
            ++cmsg ;
        }
    }
}

/******************************************/

void getPackDestSrcIdx(int dx, int dy, int dz, int fi, int i, int j, int msgType, int caseId, int &destIdx, int &srcIdx)
{
    switch (caseId) {
        case 1:
            /* Case 1 */
            srcIdx = i;
            destIdx = i;
            break;
        case 2:
            /* Case 2 */
            srcIdx = dx*dy*(dz - 1) + i;
            destIdx = i;
            break;
        case 3:
            /* Case 3 */
            srcIdx = i*dx*dy + j; 
            destIdx = i*dx + j;
            break;
        case 4:
            /* Case 4 */
            srcIdx = dx*(dy - 1) + i*dx*dy + j; 
            destIdx = i*dx + j;
            break;
        case 5:
            /* Case 5 */
            srcIdx = i*dx*dy + j*dx; 
            destIdx = i*dy + j;
            break;
        case 6:
            /* Case 6 */
            srcIdx = dx - 1 + i*dx*dy + j*dx; 
            destIdx = i*dy + j;
            break;
        case 7:
            /* Case 7 */
            srcIdx = i*dx*dy; 
            destIdx = i;
            break;
        case 8:
            /* Case 8 */
            srcIdx = i;
            destIdx = i;
            break;
        case 9:
            /* Case 9 */
            srcIdx = i*dx;
            destIdx = i;
            break;
        case 10:
            /* Case 10 */
            srcIdx = dx*dy - 1 + i*dx*dy;
            destIdx = i;
            break;
        case 11:
            /* Case 11 */
            srcIdx = dx*(dy-1) + dx*dy*(dz-1) + i; 
            destIdx = i;
            break;
        case 12:
            /* Case 12 */
            srcIdx = dx*dy*(dz-1) + dx - 1 + i*dx; 
            destIdx = i;
            break;
        case 13:
            /* Case 13 */
            srcIdx = dx*(dy-1) + i*dx*dy; 
            destIdx = i;
            break;
        case 14:
            /* Case 14 */
            srcIdx = dx*dy*(dz-1) + i; 
            destIdx = i;
            break;
        case 15:
            /* Case 15 */
            srcIdx = dx*dy*(dz-1) + i*dx; 
            destIdx = i;
            break;
        case 16:
            /* Case 16 */
            srcIdx = dx - 1 + i*dx*dy; 
            destIdx = i;
            break;
        case 17:
            /* Case 17 */
            srcIdx = dx*(dy - 1) + i; 
            destIdx = i;
            break;
        case 18:
            /* Case 18 */
            srcIdx = dx - 1 + i*dx; 
            destIdx = i;
            break;
        case 19:
            /* Case 19 */
            srcIdx = 0; 
            destIdx = fi;
            break;
        case 20:
            /* Case 20 */
            srcIdx = dx*dy*(dz - 1); 
            destIdx = fi;
            break;
        case 21:
            /* Case 21 */
            srcIdx = dx - 1; 
            destIdx = fi;
            break;
        case 22:
            /* Case 22 */
            srcIdx = dx*dy*(dz - 1) + (dx - 1); 
            destIdx = fi;
            break;
        case 23:
            /* Case 23 */
            srcIdx = dx*(dy - 1);
            destIdx = fi;
            break;
        case 24:
            /* Case 24 */
            srcIdx = dx*dy*(dz - 1) + dx*(dy - 1); 
            destIdx = fi;
            break;
        case 25:
            /* Case 25 */
            srcIdx = dx*dy - 1; 
            destIdx = fi;
            break;
        case 26:
            /* Case 26 */
            srcIdx = dx*dy*dz - 1; 
            destIdx = fi;
            break;
        default:
            assert(0 && "Invalid caseId");
    }
}

/******************************************/

void packCommon(Domain &domain, 
        int &pmsg, int &emsg, int &cmsg, 
        int xferFields, Domain_member *fieldData, 
        int n1, int n2, int msgType, int caseId, 
        int dx, int dy, int dz, 
        int destAddrIncrement) 
{
    Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
    Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
    Real_t *destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL];
    Real_t **src = (Real_t **) alloca(sizeof(Real_t *)*xferFields);
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
        src[fi] = &((domain.*fieldData[fi])(0)); 
    }

    //std::cout << "pack: msgType " << msgType << ", Case " << caseId << ". xferFields: " << xferFields << ", n1: " << n1 << ", n2: " << n2 << std::endl;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
        Domain_member src = fieldData[fi] ;
        for (Index_t i=0; i<n1; ++i) {
            for (Index_t j=0; j<n2; ++j) {
                int destIdx = -1;
                int srcIdx = -1;
                getPackDestSrcIdx(dx, dy, dz, fi, i, j, msgType, caseId, destIdx, srcIdx);
                destAddr[destIdx] = (domain.*src)(srcIdx) ;
                //std::cout << "destAddr[" << destIdx << "](&" << &destAddr[destIdx] << ")(" << destAddr[destIdx] << ") = src[" << fi << "][" << srcIdx << "](" << (domain.*src)(srcIdx) << ")" << std::endl;
            }
        }
        destAddr += destAddrIncrement ;
    }
}

/******************************************/

void CommSend(Domain& domain, int msgType,
        Index_t xferFields, Domain_member *fieldData,
        Index_t dx, Index_t dy, Index_t dz, bool doSend, bool planeOnly, int it)
{

    if (domain.numRanks() == 1)
        return ;

    //std::cout << "CommSend, msgType " << msgType << ", xferFields " << xferFields << std::endl; 

    /* post recieve buffers for all incoming messages */
    int myRank ;
    Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
    Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
    Index_t pmsg = 0 ; /* plane comm msg */
    Index_t emsg = 0 ; /* edge comm msg */
    Index_t cmsg = 0 ; /* corner comm msg */
    MPI_Datatype baseType = ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE) ;
    Real_t *destAddr ;
    bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
    /* assume communication to 6 neighbors by default */
    rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
    if (domain.rowLoc() == 0) {
        rowMin = false ;
    }
    if (domain.rowLoc() == (domain.tp()-1)) {
        rowMax = false ;
    }
    if (domain.colLoc() == 0) {
        colMin = false ;
    }
    if (domain.colLoc() == (domain.tp()-1)) {
        colMax = false ;
    }
    if (domain.planeLoc() == 0) {
        planeMin = false ;
    }
    if (domain.planeLoc() == (domain.tp()-1)) {
        planeMax = false ;
    }

    for (Index_t i=0; i<26; ++i) {
        domain.sendRequest[26*((msgType/8192)-1)+(i)] = MPI_REQUEST_NULL ;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

    /* post sends */

    if (planeMin | planeMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        int sendCount = dx * dy ;

        if (planeMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 1;
            int n1 = sendCount;
            int n2 = 1;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank - domain.tp()*domain.tp();

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 1) ;
            ++pmsg ;
        }
        if (planeMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 2;
            int n1 = sendCount;
            int n2 = 1;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank + domain.tp()*domain.tp();

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 2) ;
            ++pmsg ;
        }
    }
    if (rowMin | rowMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        int sendCount = dx * dz ;

        if (rowMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 3;
            int n1 = dz;
            int n2 = dx;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank - domain.tp();

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 3) ;
            ++pmsg ;
        }
        if (rowMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 4;
            int n1 = dz;
            int n2 = dx;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank + domain.tp();

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 4) ;
            ++pmsg ;
        }
    }
    if (colMin | colMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        int sendCount = dy * dz ;

        if (colMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 5;
            int n1 = dz;
            int n2 = dy;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank - 1;

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 5) ;
            ++pmsg ;
        }
        if (colMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm] ;
            int caseId = 6;
            int n1 = dz;
            int n2 = dy;
            int destAddrIncrement = sendCount;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int destRank = myRank + 1;

            SendCommon(destAddr, xferFields*sendCount, baseType,
                    destRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg)], 6) ;
            ++pmsg ;
        }
    }

    if (!planeOnly) {
        if (rowMin && colMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 7;
            int n1 = dz;
            int n2 = 1;
            int destAddrIncrement = dz;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp() - 1 ;

            SendCommon(destAddr, xferFields*dz, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 7) ;
            ++emsg ;
        }

        if (rowMin && planeMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 8;
            int n1 = dx;
            int n2 = 1;
            int destAddrIncrement = dx;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp()*domain.tp() - domain.tp() ;

            SendCommon(destAddr, xferFields*dx, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 8) ;
            ++emsg ;
        }

        if (colMin && planeMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 9;
            int n1 = dy;
            int n2 = 1;
            int destAddrIncrement = dy;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp()*domain.tp() - 1 ;

            SendCommon(destAddr, xferFields*dy, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 9) ;
            ++emsg ;
        }

        if (rowMax && colMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 10;
            int n1 = dz;
            int n2 = 1;
            int destAddrIncrement = dz;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp() + 1 ;

            SendCommon(destAddr, xferFields*dz, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 10) ;
            ++emsg ;
        }

        if (rowMax && planeMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 11;
            int n1 = dx;
            int n2 = 1;
            int destAddrIncrement = dx;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp()*domain.tp() + domain.tp() ;

            SendCommon(destAddr, xferFields*dx, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 11) ;
            ++emsg ;
        }

        if (colMax && planeMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 12;
            int n1 = dy;
            int n2 = 1;
            int destAddrIncrement = dy;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp()*domain.tp() + 1 ;

            SendCommon(destAddr, xferFields*dy, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 12) ;
            ++emsg ;
        }

        if (rowMax && colMin && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 13;
            int n1 = dz;
            int n2 = 1;
            int destAddrIncrement = dz;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp() - 1 ;

            SendCommon(destAddr, xferFields*dz, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 13) ;
            ++emsg ;
        }

        if (rowMin && planeMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 14;
            int n1 = dx;
            int n2 = 1;
            int destAddrIncrement = dx;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp()*domain.tp() - domain.tp() ;

            SendCommon(destAddr, xferFields*dx, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 14) ;
            ++emsg ;
        }

        if (colMin && planeMax && doSend) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 15;
            int n1 = dy;
            int n2 = 1;
            int destAddrIncrement = dy;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank + domain.tp()*domain.tp() - 1 ;

            SendCommon(destAddr, xferFields*dy, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 15) ;
            ++emsg ;
        }

        if (rowMin && colMax) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 16;
            int n1 = dz;
            int n2 = 1;
            int destAddrIncrement = dz;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp() + 1 ;

            SendCommon(destAddr, xferFields*dz, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 16) ;
            ++emsg ;
        }

        if (rowMax && planeMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 17;
            int n1 = dx;
            int n2 = 1;
            int destAddrIncrement = dx;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp()*domain.tp() + domain.tp() ;

            SendCommon(destAddr, xferFields*dx, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 17) ;
            ++emsg ;
        }

        if (colMax && planeMin) {
            destAddr = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm] ;
            int caseId = 18;
            int n1 = dy;
            int n2 = 1;
            int destAddrIncrement = dy;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            int toRank = myRank - domain.tp()*domain.tp() + 1 ;

            SendCommon(destAddr, xferFields*dy, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg)], 18) ;
            ++emsg ;
        }

        if (rowMin && colMin && planeMin) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 19;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (0, 0, 0) */
            int toRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 19) ;
            ++cmsg ;
        }
        if (rowMin && colMin && planeMax && doSend) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 20;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (0, 0, 1) */
            int toRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 20) ;
            ++cmsg ;
        }
        if (rowMin && colMax && planeMin) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 21;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (1, 0, 0) */
            int toRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 21) ;
            ++cmsg ;
        }
        if (rowMin && colMax && planeMax && doSend) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 22;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (1, 0, 1) */
            int toRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 22) ;
            ++cmsg ;
        }
        if (rowMax && colMin && planeMin) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 23;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (0, 1, 0) */
            int toRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 23) ;
            ++cmsg ;
        }
        if (rowMax && colMin && planeMax && doSend) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 24;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (0, 1, 1) */
            int toRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 24) ;
            ++cmsg ;
        }
        if (rowMax && colMax && planeMin) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 25;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (1, 1, 0) */
            int toRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 25) ;
            ++cmsg ;
        }
        if (rowMax && colMax && planeMax && doSend) {
            Real_t *comBuf = &domain.commDataSend[((msgType/8192)-1)][pmsg * maxPlaneComm +
                emsg * maxEdgeComm +
                cmsg * CACHE_COHERENCE_PAD_REAL] ;
            int caseId = 26;
            int n1 = 1;
            int n2 = 1;
            int destAddrIncrement = 0;
            packCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, destAddrIncrement); 

            /* corner at domain logical coord (1, 1, 1) */
            int toRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;

            SendCommon(comBuf, xferFields, baseType, toRank, msgType,
                    MPI_COMM_WORLD, &domain.sendRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], 26) ;
            ++cmsg ;
        }
    }
    MPI_Waitall(26, domain.sendRequest, MPI_STATUS_IGNORE) ;
}

/******************************************/

void getUnpackDestSrcIdx(int dx, int dy, int dz, int fi, int i, int j, int msgType, int caseId, int *fieldOffset, int &destIdx, int &srcIdx)
{
    if (msgType == MSG_MONOQ) {
        destIdx = fieldOffset[fi] + i; 
        srcIdx = i;
        //std::cout << "fieldOffset[" << fi << "] = " << fieldOffset[fi] << std::endl;
        //std::cout << "getUnpackDestSrcIdx for case " << caseId << " provides srcIdx " << srcIdx << ", destIdx " << destIdx << std::endl;
        return;
    }

    switch (caseId) {
        case 1:
            /* Case 1 */
            destIdx = i;
            srcIdx = i;
            break;
        case 2:
            /* Case 2 */
            destIdx = dx*dy*(dz - 1) + i;
            srcIdx = i;
            break;
        case 3:
            /* Case 3 */
            destIdx = i*dx*dy + j; 
            srcIdx = i*dx + j;
            break;
        case 4:
            /* Case 4 */
            destIdx = dx*(dy - 1) + i*dx*dy + j; 
            srcIdx = i*dx + j;
            break;
        case 5:
            /* Case 5 */
            destIdx = i*dx*dy + j*dx; 
            srcIdx = i*dy + j;
            break;
        case 6:
            /* Case 6 */
            destIdx = dx - 1 + i*dx*dy + j*dx; 
            srcIdx = i*dy + j;
            break;
        case 7:
            /* Case 7 */
            destIdx = i*dx*dy; 
            srcIdx = i;
            break;
        case 8:
            /* Case 8 */
            destIdx = i;
            srcIdx = i;
            break;
        case 9:
            /* Case 9 */
            destIdx = i*dx;
            srcIdx = i;
            break;
        case 10:
            /* Case 10 */
            destIdx = dx*dy - 1 + i*dx*dy;
            srcIdx = i;
            break;
        case 11:
            /* Case 11 */
            destIdx = dx*(dy-1) + dx*dy*(dz-1) + i; 
            srcIdx = i;
            break;
        case 12:
            /* Case 12 */
            destIdx = dx*dy*(dz-1) + dx - 1 + i*dx; 
            srcIdx = i;
            break;
        case 13:
            /* Case 13 */
            destIdx = dx*(dy-1) + i*dx*dy; 
            srcIdx = i;
            break;
        case 14:
            /* Case 14 */
            destIdx = dx*dy*(dz-1) + i; 
            srcIdx = i;
            break;
        case 15:
            /* Case 15 */
            destIdx = dx*dy*(dz-1) + i*dx; 
            srcIdx = i;
            break;
        case 16:
            /* Case 16 */
            destIdx = dx - 1 + i*dx*dy; 
            srcIdx = i;
            break;
        case 17:
            /* Case 17 */
            destIdx = dx*(dy - 1) + i; 
            srcIdx = i;
            break;
        case 18:
            /* Case 18 */
            destIdx = dx - 1 + i*dx; 
            srcIdx = i;
            break;
        case 19:
            /* Case 19 */
            destIdx = 0; 
            srcIdx = fi;
            break;
        case 20:
            /* Case 20 */
            destIdx = dx*dy*(dz - 1); 
            srcIdx = fi;
            break;
        case 21:
            /* Case 21 */
            destIdx = dx - 1; 
            srcIdx = fi;
            break;
        case 22:
            /* Case 22 */
            destIdx = dx*dy*(dz - 1) + (dx - 1); 
            srcIdx = fi;
            break;
        case 23:
            /* Case 23 */
            destIdx = dx*(dy - 1);
            srcIdx = fi;
            break;
        case 24:
            /* Case 24 */
            destIdx = dx*dy*(dz - 1) + dx*(dy - 1); 
            srcIdx = fi;
            break;
        case 25:
            /* Case 25 */
            destIdx = dx*dy - 1; 
            srcIdx = fi;
            break;
        case 26:
            /* Case 26 */
            destIdx = dx*dy*dz - 1; 
            srcIdx = fi;
            break;
        default:
            assert(0 && "Invalid caseId");
    }
}

/******************************************/

void unpackCommon(Domain &domain, 
        int &pmsg, int &emsg, int &cmsg, 
        int xferFields, Domain_member *fieldData, 
        int n1, int n2, int msgType, int caseId, 
        int dx, int dy, int dz, 
        int srcAddrIncrement, Index_t *fieldOffset = NULL) 
{
    Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
    Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
    Real_t *srcAddr = &domain.commDataRecv[((msgType/8192)-1)][pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL];
    Real_t **dest = (Real_t **) alloca(sizeof(Real_t *)*xferFields);
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
        dest[fi] = &((domain.*fieldData[fi])(0)); 
    }

    MPI_Wait(&domain.recvRequest[26*((msgType/8192)-1)+(pmsg+emsg+cmsg)], MPI_STATUS_IGNORE) ;
    //std::cerr << "Req " << 26*((msgType/8192)-1)+(pmsg+emsg+cmsg) << " completed" << std::endl;
    //std::cerr << "unpack: msgType " << msgType << ", Case " << caseId << ". xferFields: " << xferFields << ", n1: " << n1 << ", n2: " << n2 << ", srcAddr: " << srcAddr << std::endl;
    //std::cerr << "unpack: msgType " << msgType << ", Case " << caseId << ". xferFields: " << xferFields << ", n1: " << n1 << ", n2: " << n2 << std::endl;
    for (Index_t fi=0 ; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi] ;
        for (Index_t i=0; i<n1; ++i) {
            for (Index_t j=0; j<n2; ++j) {
                int destIdx = -1;
                int srcIdx = -1;
                getUnpackDestSrcIdx(dx, dy, dz, fi, i, j, msgType, caseId, fieldOffset, destIdx, srcIdx);
                //std::cout << "dest[" << fi << "][" << destIdx << "](" << (domain.*dest)(destIdx) << ") = srcAddr[" << srcIdx << "](" << srcAddr[srcIdx] << ")" << std::endl;
                //std::cerr << "dest[" << fi << "][" << destIdx << "](" << (domain.*dest)(destIdx) << ") = srcAddr[" << srcIdx << "](&" << &srcAddr[srcIdx] << ")(" << srcAddr[srcIdx] << ")" << std::endl;
                //std::cerr << "dest[" << fi << "][" << destIdx << "](" << (domain.*dest)(destIdx) << ") = srcAddr[" << srcIdx << "](" << srcAddr[srcIdx] << ")" << std::endl;
                if (msgType == MSG_COMM_SBN) {
                    (domain.*dest)(destIdx) += srcAddr[srcIdx] ;
                } else {
                    (domain.*dest)(destIdx) = srcAddr[srcIdx] ;
                }
            }
        }
        srcAddr += srcAddrIncrement ;
    }

    if (msgType == MSG_COMM_SBN || msgType == MSG_SYNC_POS_VEL) {
        if (caseId >= 0 && caseId <= 6) {
            ++pmsg;
        } else if (caseId >= 7 && caseId <= 18 ) {
            ++emsg;
        } else if (caseId >= 19 && caseId <= 26) {
            ++cmsg;
        } else {
            assert(0 && "Invalid caseId");
        }
    } else {
        //        //std::cout << "MSG_MONOQ, caseId: " << caseId << std::endl;
        for (Index_t fi=0 ; fi<xferFields; ++fi) {
            fieldOffset[fi] += srcAddrIncrement;
        }

        if (caseId >= 0 && caseId <= 6) {
            ++pmsg;
        } else {
            assert(0 && "Invalid caseId");
        }
    }
}

/******************************************/

void CommUnpackCommon(Domain& domain, int xferFields, Domain_member *fieldData, int msgType, int it) 
{
    if (domain.numRanks() == 1)
        return ;

    int myRank ;
    Index_t pmsg = 0 ; /* plane comm msg */
    Index_t emsg = 0 ; /* edge comm msg */
    Index_t cmsg = 0 ; /* corner comm msg */
    Index_t dx = msgType == MSG_MONOQ ? domain.sizeX() : domain.sizeX() + 1 ;
    Index_t dy = msgType == MSG_MONOQ ? domain.sizeY() : domain.sizeY() + 1 ;
    Index_t dz = msgType == MSG_MONOQ ? domain.sizeZ() : domain.sizeZ() + 1 ;
    bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
    bool doRecv = msgType == MSG_SYNC_POS_VEL ? false : true ;

    Index_t fieldOffset[3]; 
    fieldOffset[0] = domain.numElem() ;
    fieldOffset[1] = domain.numElem() ;
    fieldOffset[2] = domain.numElem() ;

    /* assume communication to 6 neighbors by default */
    rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
    if (domain.rowLoc() == 0) {
        rowMin = false ;
    }
    if (domain.rowLoc() == (domain.tp()-1)) {
        rowMax = false ;
    }
    if (domain.colLoc() == 0) {
        colMin = false ;
    }
    if (domain.colLoc() == (domain.tp()-1)) {
        colMax = false ;
    }
    if (domain.planeLoc() == 0) {
        planeMin = false ;
    }
    if (domain.planeLoc() == (domain.tp()-1)) {
        planeMax = false ;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

    if (planeMin | planeMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        Index_t opCount = dx * dy ;

        if (planeMin && doRecv) {
            int caseId = 1;
            int n1 = opCount;
            int n2 = 1;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
        if (planeMax) {
            int caseId = 2;
            int n1 = opCount;
            int n2 = 1;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
    }

    if (rowMin | rowMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        Index_t opCount = dx * dz ;

        if (rowMin && doRecv) {
            int caseId = 3;
            int n1 = msgType == MSG_MONOQ ? opCount : dz;
            int n2 = msgType == MSG_MONOQ ? 1 : dx;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
        if (rowMax) {
            int caseId = 4;
            int n1 = msgType == MSG_MONOQ ? opCount : dz;
            int n2 = msgType == MSG_MONOQ ? 1 : dx;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
    }

    if (colMin | colMax) {
        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        Index_t opCount = dy * dz ;

        if (colMin && doRecv) {
            int caseId = 5;
            int n1 = msgType == MSG_MONOQ ? opCount : dz;
            int n2 = msgType == MSG_MONOQ ? 1 : dy;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
        if (colMax) {
            int caseId = 6;
            int n1 = msgType == MSG_MONOQ ? opCount : dz;
            int n2 = msgType == MSG_MONOQ ? 1 : dy;
            int srcAddrIncrement = opCount;
            unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
        }
    }

    if (msgType == MSG_MONOQ) {
        return;
    }

    if (rowMin && colMin && doRecv) {
        int caseId = 7;
        int n1 = dz;
        int n2 = 1;
        int srcAddrIncrement = dz;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMin && planeMin && doRecv) {
        int caseId = 8;
        int n1 = dx;
        int n2 = 1;
        int srcAddrIncrement = dx;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (colMin && planeMin && doRecv) {
        int caseId = 9;
        int n1 = dy;
        int n2 = 1;
        int srcAddrIncrement = dy;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMax && colMax) {
        int caseId = 10;
        int n1 = dz;
        int n2 = 1;
        int srcAddrIncrement = dz;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMax && planeMax) {
        int caseId = 11;
        int n1 = dx;
        int n2 = 1;
        int srcAddrIncrement = dx;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (colMax && planeMax) {
        int caseId = 12;
        int n1 = dy;
        int n2 = 1;
        int srcAddrIncrement = dy;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMax && colMin) {
        int caseId = 13;
        int n1 = dz;
        int n2 = 1;
        int srcAddrIncrement = dz;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMin && planeMax) {
        int caseId = 14;
        int n1 = dx;
        int n2 = 1;
        int srcAddrIncrement = dx;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (colMin && planeMax) {
        int caseId = 15;
        int n1 = dy;
        int n2 = 1;
        int srcAddrIncrement = dy;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMin && colMax && doRecv) {
        int caseId = 16;
        int n1 = dz;
        int n2 = 1;
        int srcAddrIncrement = dz;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (rowMax && planeMin && doRecv) {
        int caseId = 17;
        int n1 = dx;
        int n2 = 1;
        int srcAddrIncrement = dx;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }

    if (colMax && planeMin && doRecv) {
        int caseId = 18;
        int n1 = dy;
        int n2 = 1;
        int srcAddrIncrement = dy;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }


    if (rowMin && colMin && planeMin && doRecv) {
        int caseId = 19;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMin && colMin && planeMax) {
        int caseId = 20;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMin && colMax && planeMin && doRecv) {
        int caseId = 21;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMin && colMax && planeMax) {
        int caseId = 22;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMax && colMin && planeMin && doRecv) {
        int caseId = 23;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMax && colMin && planeMax) {
        int caseId = 24;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMax && colMax && planeMin && doRecv) {
        int caseId = 25;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
    if (rowMax && colMax && planeMax) {
        int caseId = 26;
        int n1 = 1;
        int n2 = 1;
        int srcAddrIncrement = 0;
        unpackCommon(domain, pmsg, emsg, cmsg, xferFields, fieldData, n1, n2, msgType, caseId, dx, dy, dz, srcAddrIncrement, fieldOffset);
    }
}

#endif
