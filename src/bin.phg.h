#ifndef BIN_PHG_HDR
#define BIN_PHG_HDR

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define PHGRDHST_IsUsePHGHistory()		LbFgIsSet(PhgOptions, LBFlag0)		/* Will we use the PHG history file */
#define PHGRDHST_IsUseColHistory()		LbFgIsSet(PhgOptions, LBFlag1)		/* Will we use the Collimator history file */
#define PHGRDHST_IsUseDetHistory()		LbFgIsSet(PhgOptions, LBFlag2)		/* Will we use the Detector history file */

#define	PHGRDHST_NumFlags	3													/* Number of flags defined */

/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;

/* LOCAL GLOBALS */
static Boolean				phgrdhstCanceled;				/* Global cancelation flag */
static CollimatedPhotonsTy	phgrdhstColPhotons;				/* These are the successfully collimated photons */
static DetectedPhotonsTy	phgrdhstDetPhotons;				/* These are the successfully detected photons */
static char					phgrdhstErrStr[1024];			/* Error string storage */
static LbUsFourByte			phgrdhstNumToProc;				/* Number of histories to process */
static char					phgrdhstHistName[1024];			/* Name of history file */
static char					phgrdhstHistParamsName[1024];	/* Name of history parameters file */
static LbUsFourByte			phgrdhstArgIndex;
static ProdTblProdTblInfoTy	phgrdhstPrdTblInfo;				/* Info for initializing productivity table */
static PhoHFileHdrTy		phgrdhstHdrParams;				/* Input header */

/* PROTOTYPES */
Boolean 		phgrdhstInitialize(int argc, char *argv[]);
void			phgrdhstTerminate(void);
Boolean			phgbin(int argc, char *argv[]);
EventTy			readEvent(FILE *historyFile,
                                                  PHG_Decay *decayPtr,
                                                  PHG_DetectedPhoton *photonPtr);
EventTy			oldReadEvent(FILE *historyFile,
                             PHG_Decay *decayPtr,
                             PHG_DetectedPhoton *photonPtr,
                             Boolean isOldPhotons1,
                             Boolean isOldPhotons2);
Boolean			phgrdhstStandard(char *argv[]);
Boolean			phgrdhstCustom(char *argv[]);
void			phgbinProcessPhotons(PhoHFileHkTy *histHk, PHG_Decay *decayPtr,
                                     PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBlues,
                                     PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinks);



#endif
