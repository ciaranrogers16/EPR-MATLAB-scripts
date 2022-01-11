#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	27500
XMIN	10250.000000
XWID	5499.800000
*
*	Documentational Text:
*
TITL	'AB_05_CRDA340_05_50K_13000Gcentre_5500G_110secs_20dB_2G_att18dB_'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    11/10/20
TIME    16:39:19
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    1.3
A1SW    0.55
MWFQ    3.4962e+10
MWPW    0.0004647
AVGS    5
SPTP    0.004
RCAG    20
RCHM    1
B0MA    0.0002
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    27500
RCTC    0.00064
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     ezcw, 1.0

MetalIonsDef       False
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        13000.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     10250.000 G
SweepDirection     Up
SweepWidth         5500.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       34.962 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.4647 mW
PowerAtten         18 dB
QValue             800

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         5
NbScansDone        5
NbScansToDo        5
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       10

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           1
ConvTime           4.00 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               20 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             2.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           10.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadMode           False
ResCenter          62.5 mm
ResLength          20 mm
Resolution         27500
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          110.00 s
TimeConst          0.64 ms
TimeExp            False
TuneCaps           52

*
************************************************************
