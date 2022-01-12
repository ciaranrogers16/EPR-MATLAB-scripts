#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	188
XMIN	0.000000
XWID	1496.000000
*
*	Documentational Text:
*
TITL	'3PDEER_15K_9-649GHz_26-52_550MHz_Satcrash'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    04/23/21
TIME    10:57:15
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    1.992e-07
A1CT    0.322
B0VL    0.322
A1SW    1e-05
MWFQ    9.649174e+09
AVGS    0
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*
.DVC     acqStart, 1.0
.DVC     cwBridge, 1.0
AcqFineTuning      Never
AcqScanFTuning     No Scan
AcqSliceFTuning    No Slice
Power              0.0001992 mW
PowerAtten         60 dB
.DVC     endor, 1.0
EIEENDORFreq       14.9021800 MHz/3.5 kG
EIEIsotope         H1
EIERFSweepDir      Same
EIEStaticField     3480.000 G
EIEStaticRF        1.000 MHz
ENDORType          EIF
RF1Atten           60.000 dB
RF1FreqPos         1.000 MHz
RF1StartFreq       1.000 MHz
RF1SweepWidth      30.000 MHz
RF2Atten           60.000 dB
RF2FreqPos         1.000 MHz
RF2StartFreq       1.000 MHz
RF2SweepWidth      30.000 MHz
RFSrcMixing        Add
SumAtten           10.000 dB
SumAttenStart      6.000 dB
SumAttenWidth      74.000 dB
.DVC     fieldCtrl, 1.0
AllegroMode        False
CenterField        3220.00 G
Delay              0.0 s
FieldFlyback       On
FieldResol         10.00000
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3482.456 G
SweepDirection     Up
SweepWidth         0.1 G
.DVC     freqCounter, 1.0
FrequencyMon       9.649174 GHz
QMonitBridge       On
.DVC     ftBridge, 1.0
Attenuation        0.00 dB
ELDORAtt           0 dB
FrequencyA         9.102000 GHz
VideoBW            20 MHz
VideoGain          51 dB
.DVC     ftEpr, 1.0
AutoTimeOut        True
AveragesPerScan    50
ELDORFreqStart     9.150000 GHz
ELDORFreqWidth     0.000000 GHz
FTAcqModeSlct      Run from PulseSPEL
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     3P-DEER TauAvg t<0
PlsSPELGlbTxt      \
;\n\
\n\
;       Manchester Dipolar package\n\
;\n\
\n\
;       PulseSPEL general variables, definitions & conventions\n\
\n\
\n\
\n\
begin defs\n\
\n\
\n\
\n\
; Variables\n\
\n\
; Comments:\n\
\n\
\n\
\n\
; Pulse Lengths\n\
\n\
p0   = 26      ; Detection 90 pulse length\n\
\n\
p1   = 52      ; Detection 180 pulse length\n\
\n\
p2   = 0       ; Not used\n\
\n\
p3   = 20      ; "Pump" pulse length\n\
\n\
\n\
\n\
; Detection variables\n\
\n\
pg   = 100     ; Integrator Gate Width\n\
\n\
d0   = 400     ; Detection Trigger offset time 380 for 4P 400 for 3P\n\
\n\
d9   = 0       ; DAF (Delay After Flash) delay pulse sequence until\n\
\n\
                    ; after external event (eg. laser pulse)\n\
\n\
\n\
\n\
; Dipolar Delay Times\n\
\n\
d1   = 400     ; tau1\n\
\n\
d2   = 1500    ; tau2\n\
\n\
d3   = 300     ; initial T in 4P DEER d1 - d3  = time before zero time\n\
\n\
d4   = 100     ; time before zero time in 3P DEER\n\
\n\
\n\
\n\
; Dipolar Time Steps\n\
\n\
d30  = 8       ; 1st time scale increment\n\
\n\
d31  = 56      ; time step for tau averaging\n\
\n\
p    = 8       ; number of tau1s to average\n\
q    = 8       ; number of tau2s to average\n\
\n\
; Tuning/Nutation Delay Times\n\
\n\
d11  = 400     ; Delay in Hahn Echo detection\n\
\n\
d12  = 5000    ; Delay after inversion pulse\n\
\n\
; Nutation Time Step\n\
\n\
d29  = 2       ; Length increment of `inversion' pulse\n\
\n\
\n\
\n\
; Experiment Parameters\n\
\n\
a    = 64      ; (RESERVED) number of transient averages (TRANS)\n\
\n\
                ; if greater than 1024, a must be a multiple of 1024\n\
\n\
h    = 50       ; (CONVENTION) number of shots per point (INTG and SPT)\n\
\n\
n    = 400     ; (CONVENTION) number of sweeps to accumulate\n\
\n\
\n\
\n\
srt  = 3000 * srtu ; (RESERVED) SRT - Shot Repetition Time (srtu=1.02 us)\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
; Scared to remove\n\
\n\
b    = 1       ;\n\
\n\
c    = 1       ;\n\
\n\
r    = 1       ;\n\
\n\
s    = 1       ;\n\
\n\
t    = 1       ;\n\
\n\
\n\
\n\
w    = 800     ; (RESERVED) Split point of sweeps (ns)\n\
\n\
\n\
\n\
dx   = 0       ; (RESERVED) 1st time scale starting value\n\
\n\
dy   = 0       ; (RESERVED) 2nd time scale starting value\n\
\n\
\n\
\n\
end defs\n\
\n\
\n\
\n\
; Note:\n\
\n\
;       Variables Marked 'RESERVED' = variable may be used for other times,\n\
\n\
;               however only this variable can change\n\
\n\
;               the indicated parameter.\n\
\n\
;               eg. only 'a' will change the number of\n\
\n\
;               transients averaged in SpecJet\n\
\n\
\n
PlsSPELLISTSlct    DEER-2 dMPFU/pELDOR
PlsSPELPhPrgEx     Normal
PlsSPELPrg         PulseSPEL/AliceDEERMIB.exp
PlsSPELPrgTxt      \
;\n\
; 4P-DEER and associated tuning routines\n\
;\n\
begin defs              ;       dimension of data array [sx,sy]\n\
 dim    s[512]          ;       for MPFU HE detection 90 pulse +/-<x>\n\
 dim1   s[512]          ;       for ELDOR pulse nutation\n\
 dim2   s[1024]         ;       for 4P-DEER phase-cycled transient\n\
 dim3   s[188]          ;       for 4P-DEER with 2D tau averaging  ca. (d1+d2-d3-80)/d30 this assumes your pump pulse is <80 ns\n\
 dim4   s[188]		;       for 3P-DEER with tau averaging and t<0  ca. (d4+d2-80)/d30 this assumes your pump pulse is <80 ns\n\
end defs\n\
\n\
\n\
begin lists "+<x>(pi/2)"\n\
 ph1 +<x>\n\
 ph2 +<x>\n\
 ph3 *\n\
 ph4 *\n\
 asg1 +a\n\
 bsg1 +b\n\
end lists \n\
\n\
begin lists2 "-<x>(pi/2)"\n\
 ph1 -<x>  \n\
 ph2 +<x>\n\
 ph3 *\n\
 ph4 *\n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists2\n\
\n\
\n\
begin lists3 "nutELDOR d+<x>"\n\
 ph1 +<x>   \n\
 ph2 +<x>\n\
 ph3 ELDOR\n\
 ph4 *\n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists3\n\
\n\
\n\
begin lists4 "DEER-2 dMPFU/pELDOR"\n\
 ph1 +<x> -<x>\n\
 ph2 +<x>\n\
 ph3 ELDOR\n\
 ph4 +<x>             \n\
 asg1 +a -a                  \n\
 bsg1 +b -b\n\
end lists4\n\
\n\
\n\
begin lists5 "DEER-4 dMPFU/pELDOR"\n\
 ph1 +<x> -<x>\n\
 ph2 +<x> +<x> -<x> -<x>\n\
 ph3 ELDOR\n\
 ph4 +<x> +<x> -<x> -<x>\n\
 asg1 +a -a +a -a                  \n\
 bsg1 +b -b +b -b\n\
end lists5\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
\n\
; Experiment Definitions\n\
\n\
;\n\
; MPFU HE detection 90 pulse +/-<x>\n\
;\n\
begin exp "MPFU 90 pulse +/-<x> : +<x> 1st" [TRANS QUAD]\n\
; QUAD detection with TRANSient visulalization\n\
\n\
p5 = 2*p0                       ; pi pulse length\n\
    shot i = 1 to h             ; accumulation loop\n\
        d9                      ; DAF\n\
        p0 [ph1]                ; pi/2\n\
        d11                     ; tau1\n\
        p5 [ph2]                ; pi\n\
        d11                     ; tau1\n\
        d0                      ; acquisition delay\n\
        dig [sg1]               ; acquisition\n\
    next i                      ; end of accumulation loop\n\
\n\
end exp\n\
\n\
\n\
\n\
;\n\
; MPFU HE detection 180 pulse +/-<x>\n\
;\n\
begin exp1 "ELDOR nutation" [TRANS QUAD]\n\
; QUAD detection with TRANSient visualization\n\
p5 = p1/2                       ; pi/2 pulse length\n\
    shot i = 1 to h             ; accumulation loop\n\
        p3 [ph3]\n\
        d12                     ; DAF\n\
        p5 [ph1]                ; pi/2\n\
        d11                     ; tau1\n\
        p1 [ph2]                ; pi\n\
        d11                     ; tau1\n\
        d0                      ; acquisition delay\n\
        dig [sg1]               ; acquisition\n\
    next i                      ; end of accumulation loop\n\
\n\
end exp1\n\
\n\
\n\
;\n\
; 4P-DEER Transient\n\
;\n\
begin exp2 "4P-DEER Transient"  [TRANS QUAD]\n\
; QUAD detection of Transient\n\
\n\
for k = 1 to n                  ; averaging loop\n\
  totscans(v)                   ; output of total number of scans\n\
  shot i = 1 to a               ; accumulation loop\n\
        d9                      ; DAF (Delay After Flash)\n\
        p0 [ph1]                ; 1st pulse\n\
        d1                      ; tau1\n\
        p1 [ph2]                ; 2nd pulse\n\
        d1                      ; delay to "pump"\n\
        p3 [ph3]                ; "pump"\n\
        d2                      ; delay after "pump"\n\
        p1 [ph4]                ; 4th pulse\n\
        d2                      ; tau2\n\
        d0                      ; constant acquisition delay\n\
        dig [sg1]               ; signal acquisition\n\
  next i                        ; end of accumulation loop\n\
  scansdone(k)          ; output number of scans done\n\
next k\n\
\n\
end exp2\n\
\n\
\n\
\n\
;\n\
; 4P-DEER with 2D Tau averaging onboard\n\
;\n\
begin exp3 "4P-DEER 2DTauAvg"  [INTG QUAD]\n\
; QUAD detection with Integration\n\
\n\
d21 = d1                        ; store tau1\n\
d22 = d2                        ; store tau2\n\
d23 = d3                        ; store initial T\n\
\n\
v = p * q                       ; number of tau averaging steps\n\
v = v * n                       ; total number of scans to run\n\
\n\
for k = 1 to n                  ; averaging loop\n\
  d1 = d21                      ; reset tau1\n\
  d3 = d23                      ; reset delay before "pump"\n\
  v = p * q                     ; number of tau averaging steps\n\
  v = v * n                     ; total number of scans to run\n\
  totscans(v)                   ; output of total number of scans\n\
  for j = 1 to p                ; loop over tau1 steps\n\
    d2 = d22                    ; reset tau2\n\
    for l = 1 to q              ; loop over tau2 steps\n\
      d8 = d3                   ; reset delay before "pump"\n\
      d7 = d1 + d2              ; total delay time available\n\
      d7 = d7 - d8              ; reset delay after "pump"\n\
      dx = d8 - d1              ; reset x-axis\n\
      sweep x = 1 to sx         ; T sweep loop\n\
       shot i = 1 to h          ; accumulation loop\n\
        d9                      ; DAF (Delay After Flash)\n\
        p0 [ph1]                ; 1st pulse\n\
        d1                      ; tau1\n\
        p1 [ph2]                ; 2nd pulse\n\
        d8                      ; delay to "pump"\n\
        p3 [ph3]                ; "pump"\n\
        d7                      ; delay after "pump"\n\
        p1 [ph4]                ; 4th pulse\n\
        d2                      ; tau2\n\
        d0                      ; constant acquisition delay\n\
        acq [sg1]               ; signal acquisition\n\
       next i                   ; end of accumulation loop\n\
       d8 = d8 + d30            ; increment time before "pump"\n\
       d7 = d7 - d30            ; decrement time after "pump"\n\
       dx = dx + d30            ; increment x axis\n\
      next x                    ; end of T sweep loop\n\
      d2 = d2 + d31             ; increment tau2\n\
      u = u + 1                 ; add to total scans completed\n\
      scansdone(u)              ; output number of scans done\n\
    next l                      ; next tau2-loop\n\
    d1 = d1 + d31               ; increment tau1\n\
    d3 = d3 + d31               ; increment initial delay before "pump"\n\
  next j                        ; next tau1-loop\n\
next k\n\
\n\
end exp3\n\
\n\
\n\
;\n\
; 4P-DEER with 2D Tau averaging onboard\n\
;\n\
begin exp4 "3P-DEER TauAvg t<0"  [INTG QUAD]\n\
; QUAD detection with Integration\n\
\n\
l = d4/d30    			; number of points for pump pulse before 1st det. pulse\n\
b = l+1				; point of aligned pump and det. pulses\n\
u = 0\n\
\n\
for k = 1 to n                  ; averaging loop\n\
\n\
v = p * n                       ; total number of scans to run using tau averaging\n\
totscans(v)			; output number of scans\n\
\n\
d6 = d2                         ; store tau\n\
\n\
for j = 1 to p 			; tau averaging loop\n\
; record DEER trace with pump pulse before detection\n\
d7 = 0				; set initial time before pump (delay 1) note: d7 + d8 = d4 (constant)\n\
d8 = d4				; set initial time adter pump before first det pulse (delay 2)\n\
dx = 0				; define x axis\n\
\n\
sweep x = 1 to l\n\
   shot i = 1 to h\n\
	d9                      ; DAF (Delay After Flash)\n\
	d7			; varialbe pump delay 1\n\
        p3 [ph3]                ; "pump"\n\
        d8			; varialbe pump delay 2\n\
        p0 [ph1]                ; 1st det pulse\n\
        d6                      ; tau\n\
        p1 [ph2]                ; 2nd det pulse\n\
        d6                      ; tau\n\
        d0                      ; constant acquisition delay\n\
        acq [sg1]               ; signal acquisition\n\
     next i                  	; end of accumulation loop\n\
     d7 = d7 + d30           	; increment time before "pump"\n\
     d8 = d8 - d30           	; decrement time after "pump"\n\
     dx = dx + d30          	; increment x axis\n\
next x                    	; end of T sweep loop\n\
\n\
\n\
; record DEER trace with pump pulse after first detection pulse\n\
d7 = 0				; set initial time before pump\n\
d8 = d6				; set initial time after pump before pi det.\n\
\n\
sweep x = b to sx\n\
   shot i = 1 to h\n\
	d9                      ; DAF (Delay After Flash)\n\
        p0 [ph1]                ; 1st det pulse\n\
	d7			; variable pump delay\n\
        p3 [ph3]                ; "pump"\n\
        d8                      ; variable pump delay\n\
        p1 [ph2]                ; 2nd det pulse\n\
        d6                      ; tau\n\
        d0                      ; constant acquisition delay\n\
        acq [sg1]               ; signal acquisition\n\
     next i                     ; end of accumulation loop\n\
     d7 = d7 + d30              ; increment time before "pump"\n\
     d8 = d8 - d30              ; decrement time after "pump"\n\
     dx = dx + d30              ; increment x axis\n\
next x                          ; end of T sweep loop \n\
d6 = d6 + d31			; tau averaging increment \n\
u=u+1\n\
scansdone(u)   \n\
next j\n\
\n\
next k\n\
\n\
end exp4\n\
\n
Psd1               0,0,276,1776,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,30,56,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              0,268,268,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,28,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1836,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1736,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               0,0,3400,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               0,264,264,1764,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,28,54,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               0,276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               0,264,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        3060.00 us
ShotsPLoop         50
SptProgress        21 %
StochMode          False
SweepsPExp         3200
TriggerTimeOut     41 s
XAxisQuant         Time
XSpecRes           188
YAxisQuant         Magnetic Field
YSpecRes           1
.DVC     recorder, 1.0
BaselineCorr       Off
NbScansAcc         0
NbScansDone        0
NbScansToDo        3200
ReplaceMode        Off
.DVC     routeTrg, 1.0
.DVC     sigChanSmall, 1.0
.DVC     transRec, 1.0
*
************************************************************
