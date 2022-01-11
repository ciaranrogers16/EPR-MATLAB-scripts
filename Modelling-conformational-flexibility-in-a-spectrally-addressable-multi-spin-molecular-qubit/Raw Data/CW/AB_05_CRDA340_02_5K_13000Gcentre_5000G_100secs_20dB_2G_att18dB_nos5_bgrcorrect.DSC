#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	MAN
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
XPTS	25000
XMIN	10500.000000
XWID	4999.800000
*
*	Documentational Text:
*
TITL	'AB_05_CRDA340_02_5K_13000Gcentre_5000G_100secs_20dB_2G_att18dB_n'
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
TIME    15:42:46
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
A1SW    0.5
MWFQ    3.4937e+10
MWPW    0.0006963
AVGS    5
SPTP    0.004
RCAG    20
RCHM    1
B0MA    0.0002
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    25000
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
StaticFieldMon     10500.000 G
SweepDirection     Up
SweepWidth         5000.0 G
WidthTM            200.0 G
.DVC     fieldSweep, 1.0
.DVC     freqCounter, 1.0
FrequencyMon       34.937 GHz
QMonitBridge       On
.DVC     mwBridge, 1.0
AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.6963 mW
PowerAtten         18 dB
QValue             600
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
Resolution         25000
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          100.00 s
TimeConst          0.64 ms
TimeExp            False
TuneCaps           52
*
************************************************************
*
#MHL	1.0 * MANIPULATION HISTORY LAYER by BRUKER
*
*
*       Input Dataset:
*
SOURCE_PRIM
'AB_05_CRDA340_02_5K_13000Gcentre_5000G_100secs_20dB_2G_att18dB_n'
END_SOURCE_PRIM
*
*       Input Dataset:
*
SOURCE_SCND
  *
  SOURCE
  'AB_05_CRDA340_02_5K_13000Gcentre_5000G_100secs_20dB_2G_att18dB_n'
  END_SOURCE
  SELECT  'qualiRegions'
  'qualiRegions'
  END_SELECT
  PROCESS 'prLinRegr'
  REPORT
                  a = -0.0000831 +- 0.0000001
                  b = -3.537     +- 0.002
  reduced chi-square = 0.0002607
  END_REPORT
END_SOURCE_SCND
*
*       Selection:
*
SELECT  'qualiRegions'
'qualiRegions'
END_SELECT
*
*       Process:
*
PROCESS 'prDiff'
*
*       Parameter:
*
PAR_VAL Gain(Sec.) =  1.000e+00
PAR_VAL x-Shift(Sec.) =  0.000e+00
PAR_VAL x-Stretch(Sec.) =  1.000e+00
MDATE   11/12/20 16:19:37
*
************************************************************
