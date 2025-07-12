!--------+---------+---------+---------+---------+---------+---------+--
PROGRAM p17                                                     
use modu
!                                                                       
! ==================================================================================
!     DICHIARAZIONE VARIABILI                      
! ==================================================================================
! 
REAL,ALLOCATABLE::  QF(:,:,:),QC(:,:,:),DS(:,:,:),CR(:,:,:),DL(:,:),DD(:,:),        &
		    CP(:,:),DDUBJ(:),SIGMA(:,:),DUB(:,:),DUBW(:,:),A(:,:),RHS(:),   &
		    A1(:,:),RHS0(:),QFW(:,:,:),DQFW(:,:,:),dubprec(:,:),PHISIG(:,:),&
		    PHISIGPREC(:,:),UBC(:),WBC(:)
REAL::		    ZERO,T1
INTEGER,ALLOCATABLE:: IP(:) 
!                                                                       
COMMON/NO1/ CROOT,CTIP,XTIP,ZTIP,B,S,AR,PAY 
COMMON/NO2/ RO,Q_INF,UT 
COMMON/NO3/ FOURPAY 
!
CHARACTER(LEN=3)::  STRIB,STRJB,STRA,STRB,STRV,STRKH
! ----------------------------------------------------------------------------------
! ==================================================================================
!     LETTURA INPUT E ALLOCAZIONE
! ==================================================================================
!
open(1, file = 'input.dat')
    read(1,*)	IB
    read(1,*)	JB
    read(1,*)	NSTEPS
    read(1,*)	B
    read(1,*)	Q_INF
    read(1,*)	ALPHA1
    read(1,*)	TEMP
    read(1,*)	AA   
    read(1,*)	OMA  
    read(1,*)	SZA  
    read(1,*)	OMH  
    read(1,*)	FASE     
    read(1,*)	THIC 
close(1)
!
NDIM = IB * JB
IB1  = IB + 1
JB1  = JB + 1
IB2  = IB + 2
JM   = JB/2 + 1
IJB = IB * JB
!
ALLOCATE(   QF(IB2,JB1,3),QC(IB1,JB,3),DS(IB1,JB,10),CR(IB1,JB,12),DL(IB,JB),	   &
	    DD(IB,JB),CP(IB,JB),DDUBJ(JB),SIGMA(IB,JB),DUB(IB,JB),DUBW(NSTEPS,JB), &
	    A(NDIM,NDIM),RHS(NDIM),A1(IB,JB),RHS0(NDIM),QFW(NSTEPS+1,JB+1,3),      &
	    DQFW(NSTEPS,JB1,3),IP(NDIM),dubprec(ib,jb),PHISIG(IB,JB),		   &
	    PHISIGPREC(IB,JB),UBC(IB),WBC(IB)	    )
! ----------------------------------------------------------------------------------                                                                   
! ==================================================================================
!     CREAZIONE OUTPUT E INIZIALIZZAZIONE
! ==================================================================================
!
! outfile simply converts the first variable to a string character
!
CALL OUTFILE(ib,strib)
CALL OUTFILE(jb,strjb)
CALL OUTFILE(nint(b),strb)
CALL OUTFILE(nint(alpha1),stra)
CALL OUTFILE(nint(q_inf),strv)
!
open(6,	file  =	"./out/output.csv") 
open(7,	file  =	"./out/cl_vs_t.csv") 
open(9, file  = "./out/wakepoints.csv")
open(17,file  = "./out/cp_mid.csv")
open(19, file = "./out/wake_3d.csv")
open(66,file  =	"./out/wake_gnd.csv") 
!                                                                       
PAY=acos(-1.) 
FOURPAY = 4.*PAY 
RO	= 1.0 
ALFA0	= ALPHA1*PAY/180.
alfamp	= aa*pay/180. 
CROOT	= 1.0                                        
CTIP	= 1.0 
XTIP	= 0.0 !AFT SWEEP OF TIP 
ZTIP	= 0.0 
S = 0.5*B*(CROOT + CTIP)
AR = B**2/S
!                                                  
UT    = Q_INF 
CM    = (CROOT+CTIP)/2. 
DX    = CM/(IB/2) 
DLW   = 0.05 
!DXW = DLW*CSA 
!DZW = DLW*SNA 
DT  = DLW/Q_INF/temp
!                                                                       
!dt = 0.004
!     ORIGINE S.R. GROUND                                               
X0 = 0. 
Y0 = 0. 
Z0 = 0.
! 
WRITE(6,*) "N,CL,FL,CM,CD"
! WRITE(6,101) 
!   101 FORMAT(1H ,/,20X,'INTERNAL POTENTIAL BASED PANEL CODE',/,15X,46('-'))         
! WRITE(6,102) ALPHA1,B,CROOT,S,AR,Q_INF,IB,JB 
!   102 FORMAT(1H ,/,10X,'ALFA:',F10.2,8X,'B   :',F10.2,8X,'C   :',F13.2,/,10X,&
!      'S   :',F10.2,8X,'AR  :',F10.2,8X,'V(INF) :',F10.2,/,10X,&
!      'IB  :',I10,8X,'JB  :',I10,8X,/)                                  
! !     WRITE(6,111)                                                      
! write(6,*)'DX=',DX,'DT=',DT,'Q_INF*DT=',Q_INF*DT 
! ONE=1. 
! write(6,*)'PRECIS=',EPSILON(ONE) 
! WRITE(6,113) 
!WRITE(8,113) 
113 FORMAT(1H ,73('='))
!
DO NN=1,NSTEPS 
    DO J=1,JB 
        DUBW(NN,J)=0.0 
    ENDDO                                                          
ENDDO
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
T = -DT
DO N = 1,NSTEPS
! ----------------------------------------------------------------------------------
! ==================================================================================
!	CINEMATICA
! ==================================================================================
    T = T + DT
    ALFA    = ALFA0 + ALFAMP*COS(OMA*T)
    SZ      = SZA*COS(OMH*T + FASE*PAY/180.) ! DEL CORPO
    ALFADOT = -OMA*ALFAMP*SIN(OMA*T)
    SZDOT   = -OMH*SZA*SIN(OMH*T + FASE*PAY/180.)
!
    CSA = COS(ALFA)
    SNA = SIN(ALFA)
! ---------------------------------------------------------------------------------- 
! ================================================================================== 
!	WING GEOMETRY                                                     
! ================================================================================== 
!
    CALL GRID(IB,JB,SZ,QF,QC,DS,CR,DLW,CSA,SNA,thic) 
!    
    IF (N.GT.1) THEN
	PHISIGPREC = PHISIG
    ELSE
	PHISIGPREC = 0.
    END IF
!
    IM = IB2/2
    DO I = 1,IB
	DO J = 1,JB                                                     
	    UBC(I) = -UT + ALFADOT*(QC(I,J,3)-QF(IM,J,3))
	    WBC(I) = SZDOT -ALFADOT*(QC(I,J,1)-QF(IM,J,1))
	    SIGMA(I,J) = DS(I,J,7)*UBC(I) + DS(I,J,9)*WBC(I)  !n_ext \cdot v_BODY
	END DO
    END DO
!                                                                      

open(22,file= 'qf.txt')
	do i = 1,ib2
	    write(22,*) qf(i,1,1),qf(i,1,2),qf(i,1,3)
	end do
	write(22,*) '#'
!     =============                                                     
!     PROGRAM START                                                     
!     =============                                                     
!
!     INFLUENCE COEFFICIENTS CALCULATION                                
!                                                                       
               !COLLOCATION POINT COUNTER                               
    K=0 
    DO I=1,IB 
    DO J=1,JB 
	K=K+1 
	CALL WING(IB,JB,QC,CR,I,J,DDUBJ,DS,A1,SIGMA,RH) 
	L=0 
	DO I1=1,IB 
        DO J1=1,JB 
        !INFLUENCING PANEL COUNTER                 
	    L=L+1 
!		A(K,L) - IS THE VELOCITY POTENTIAL ...                    
!		DUE TO A UNIT DOUBLET                                     
		A(K,L)=A1(I1,J1) 
	ENDDO 
	ENDDO 

!                                                                       
!       CALCULATE RHS                                                 
	RHS0(K)=-RH 
	RHS(K)=RHS0(K) 
    ENDDO 
    ENDDO                                                                        

!   ===========================                                     
!   VORTEX WAKE SHEDDING POINTS                                     
!   ===========================                                     
    DO J=1,JB1 
        QFW(N,J,1)=QF(IB2,J,1) 
        QFW(N,J,2)=QF(IB2,J,2) 
        QFW(N,J,3)=QF(IB2,J,3) 
!
        QFW(N+1,J,1)=QF(IB1,J,1) 
        QFW(N+1,J,2)=QF(IB1,J,2) 
        QFW(N+1,J,3)=QF(IB1,J,3) 
    ENDDO 
!                                                                       
!   ========================                                        
!   AERODYNAMIC CALCULATIONS                                        
!   ========================                                        
!                                                                       
!       INFLUENCE COEFFICIENTS CALCULATION da scia su ALA               
!                                                                       
                 !COLLOCATION POINT COUNTER                             
    K=0 
    DO I=1,IB 
    DO J=1,JB 
	K=K+1 
!       CALCULATE WAKE INFLUENCE (P16_5b.f riga 201)                
        CALL POT_WAKE(IB,JB,NSTEPS,QFW,QC,DUBW,I,J,N,RH) 
!                                                                       
        RHS(K)=RHS0(k) - RH 
    ENDDO 
    ENDDO 

!                                                                       
!   SOLUZIONE                           
    CALL DECOMP(NDIM,NDIM,A,IP)
!   SALVATAGGIO DOPPIETTE ISTANTE PRECEDENTE
    if (n.gt.1) then
	do i = 1,ib
	do j = 1,jb
	    dubprec(i,j) = dub(i,j)
	end do 
	end do
    else
	dubprec = 0.
    end if
    CALL SOLVER(NDIM,NDIM,A,RHS,IP) 
!
!   WING DOUBLET LATTICE LISTING                                    
    K=0 
    DO I=1,IB 
	DO J=1,JB 
	    K=K+1 
            DUB(I,J)=RHS(K) 
        ENDDO 
    999       FORMAT(I3,1(1x,F9.6),A25) 
    ENDDO 
!
!                                                                       
!       WAKE SHEDDING                                                   
!                                                                       
    DO J=1,JB 
	DUBW(N,J) = DUB(IB,J)-DUB(1,J) 
    ENDDO 
!                                                              
    ! calcolo di PHI per dphi/dt
    DO I = 1,IB
    DO J = 1,JB
	PHISIG(I,J) = - DUB(I,J)
    END DO
    END DO
    CALL FORCES(IB,JB,QF,QC,DS,DUB,DUBPREC,PHISIG,PHISIGPREC,dt,UBC,WBC,N,DL,DD,CP,CL,CD,CM,FL,QI,QM,QL) 
!        
!   Write CL versus time results
!
    write(7,*)T,omh*T/(2*PAY),OMA*T/(2*PAY),OMH*T,OMA*T,CL,sz,' t OMH/2PI  OMA/2PI OMH OMA	CL'
!                                                                       
!       ===========================                                     
!       WAKE ROLLUP CALCULATION see P16_5b.f 266-291                    
!       ===========================                                     
!                                                                       
    IW=1 
!       IF(N.GT.1) THEN                                                 
!         IF(N.GE.NW) IW=N-NW+1                                         
!         NW IS THE NUMBER OF (TIMEWISE) DEFORMING WAKE ELEMENTS.       
!         deforma solo gli elementi della scia da N-NW+1 a N-1          
                                                                        
                     !-1 ! influenced wake panel                        
    DO IT=IW,N 
	!write(23,*) 
        !write(23,*)N,IT,'------- N,IT ------' 
                       ! IT,J è il generico nodo della scia QFW        
        DO J=1,JB1 
	    X = QFW(IT,J,1) 
            Y = QFW(IT,J,2) 
            Z = QFW(IT,J,3) 
! calcolo vel su scia dovuta a sorg e doppiette su corpo                
            CALL VELOCE2(IB,JB,QC,CR,DUB,SIGMA,DS,X,Y,Z,U,V,W) 
! 123   FORMAT(3(F12.7,1x),A30)                                         
! calcolo vel su scia dovuta a doppiette scia                           
            CALL VEL_WAKE(JB,NSTEPS,X,Y,Z,QFW,DUBW,N,UW,VW,WW) 
!                                                                       
            U = U + UW 
            V = V + VW 
            W = W + WW 
! calcola spostamento nodi scia nel SR corpo                            
            DQFW(IT,J,1)=(U+UT)*DT 
            DQFW(IT,J,2)=(V   )*DT 
            DQFW(IT,J,3)=(W)*DT 
        ENDDO 
    ENDDO 
!                                                                       
!                                                                       
! aggiorna posizione nodi scia al tempo nuovo (sempre nel SR corpo)     
    DO IT=1,N 
	DO J=1,JB1 
	    QFW(IT,J,1)=QFW(IT,J,1)+DQFW(IT,J,1) 
            QFW(IT,J,2)=QFW(IT,J,2)+DQFW(IT,J,2) 
            QFW(IT,J,3)=QFW(IT,J,3)+DQFW(IT,J,3) 
        ENDDO
    !    
    ! Write Wake points coordinates    
	!
    write(9,*) QFW(IT,JM,1),QFW(IT,JM,2),QFW(IT,JM,3) 
    ENDDO
    write(9,*)'#' 
!!!!
! -------------------------------------------------------------------
! ===================================================================
!	STAMPA OUTPUT
! ===================================================================
! -------------------------------------------------------------------  

    !DO IT=1,N+1 
	!write(25,133)N,IT,QFW(IT,JM,1),QFW(IT,JM,3),JM,' N IT QFW(IT,JM,1/3) JM' 
    !ENDDO 
    133   FORMAT(2(I2,1x),2(F12.7,1x),I3,A35) 
    !WRITE(25,*) 
!                                                                       
!     OUTPUT                                                            
!                    
!   write results in main output file
!                                                   
    
    WRITE(6,104) N,CL,FL,CM,CD 
    104   FORMAT(1H ,I3,",",F10.4,",",F10.4,",",F10.4,",",F10.4) 
    IF(N.EQ.1) THEN 
!       scrive filamenti streamwise rilasciati al passo 1               
	!j=1 
        !write(79,*)'# filamenti streamwise rilasciati al passo 1' 
        !write(79,*)J,QFW(N,J,2),-DUBW(N,J),'J,QFW(N,J,2),circ' 
        DO J=2,JB 
	    circ=DUBW(N,J-1)-DUBW(N,J) 
            !write(79,*)j,QFW(N,J,2),circ,'J,QFW(N,J,2),circ' 
        ENDDO 
!       scrive circolazione della riga 1 spanwise della scia            
        !j=JB1 
        !write(80,*)'# circolazione della riga 1 spanwise della scia' 
        !write(79,*)J,QFW(N,J,2), DUBW(N,JB),'J,QFW(N,J,2),circ' 
        DO J=1,JB 
	    yy=0.5*(QFW(N,J,2)+QFW(N,J+1,2)) 
            !write(80,*)j,yy,DUBW(N,J),'j,yy,DUBW(N,J)' 
        ENDDO 
    ENDIF 
            !N  AA                                                      
ENDDO ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                     
JM = JB/2
do i = 1,ib
!
! writes pressure coefficient chordwise at wing center
!
	write(17,*) QC(I,JM,1),CP(I,JM)
end do
DO I = 1,IB1
!
!writes coordinates of the wing
!
	write(19,*) (QF(I,J,1),j=1,jb1)
END DO
DO I = 1,IB1
	write(19,*) (QF(I,J,2),j=1,jb1)
END DO
DO I = 1,IB1
	write(19,*) (QF(I,J,3),j=1,jb1)
END DO
DO IT=1,NSTEPS 
!
!writes coordinates of the wakepoints
!
	write(19,*) (QFW(IT,J,1),j=1,jb1)
ENDDO
DO IT=1,NSTEPS 
	write(19,*) (QFW(IT,J,2),j=1,jb1)
ENDDO
DO IT=1,NSTEPS 
	write(19,*) (QFW(IT,J,3),j=1,jb1)
ENDDO
! 
110 FORMAT(/,1H ,2X,'I   J',7X,'XC',7X,'CP',8X,'DL',8X,'DD',7X,'DUB',6X,'SIGMA',/,68('=')) 
!       --------                                                        
!    
!
! Writes wing points and wake points in fixed frame of reference for wing animation
!
!                         
WRITE(66,*)'# coord del corpo e della scia (J=JB/2+1) nel S.R. ground'    
DO I=1,IB2 
    WRITE(66,123) QF(i,1,1),QF(I,1,2),QF(I,1,3),' QF(I,1,1-3) post' 
    123   FORMAT(3(F12.7,1x),A30)                                        
ENDDO 
WRITE(66,*)'#' 
!       sezione centrale scia                                           
DO IT=NSTEPS,1,-1 
    QF1 = QFW(IT,JB/2 +1,1)
    QF3 = QFW(IT,JB/2 +1,3) 
    WRITE(66,125) QF1,QFW(IT,JB/2 +1,2),QF3,IT,' QFW(IT,JM,x-y-z) '         
ENDDO 
125   FORMAT(3(F12.7,1x),i4,A30) 
WRITE(66,*) 
!                                                                       
      STOP 
END PROGRAM                                           
!                                                                       
                                           
