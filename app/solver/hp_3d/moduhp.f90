module modu
contains
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE BODY(IB,JB,QF,thic) 
!     DEFINISCE LA PRIMA SEZIONE DEL CORPO CILINDRICO                   
      DIMENSION QF(IB+2,JB+1,3) 
!     a1...a5 dati per profilo NACA 00t                                 
      data a1,a2,a3,a4,a5/1.4845,-.63,-1.7685,1.4215,-0.5075/ 
      a2=-(a1+a3+a4+a5) 
!                                                                       
!     i pannelli sono numerati partendo dal LE dorso e                  
!     procedendo in senso orario                                        
      pi=acos(-1.0) 
      nm=IB/2 
               !NACA 0012                                               
      t=thic 
      dth=pi/float(nm) 
      DO i=1,nm+1 
        the=(i-1.)*dth 
                              !full cos. spacing                        
        z=(-cos(the)+1.)*.5 
!       z=(i-1.)/float(nm)    !dx = cost                                
        yy=t*(a1*z**.5+a2*z+a3*z*z+a4*z**3+a5*z**4) 
        yy=t*(a1*sqrt(z)+z*(a2+z*(a3+z*(a4+a5*z)))) 
                         !da nm+1 a n+1 dorso                           
        QF(i+nm,1,1)=z 
        QF(i+nm,1,3)=yy 
                         !da nm+1 a 1 ventre                            
        QF(nm+2-i,1,1)=z 
        QF(nm+2-i,1,3)=-yy 
      ENDDO 
      QF(IB+1,1,1)=QF(1,1,1) 
      QF(IB+1,1,3)=QF(1,1,3) 
!     do i=1,IB+1                                                       
!       write(*,10) QF(I,1,1),QF(I,1,3)                                 
!     ENDDO                                                             
!  10 FORMAT(1x,F7.4,2x,F7.4)                                      
      RETURN 
      END                                           
!-----------------------------------------------------------------------
      SUBROUTINE GRID(IB,JB,sz,QF,QC,DS,CR,DXW,CSA,SNA,thic) 
!     CALCULATE PANEL CORNERPOINTS;  QF(I,J,(X,Y,Z))                    
!               WING COLLOCATION POINTS QC(I,J,(X,Y,Z))                 
!               SOURCE STRENGTH SIGMA(I,J)                              
!               TANGENT, NORMAL and AREA on EACH PANEL DS(I,J,1-10)     
!               PANEL CORNERPOINTS IN WING R.S. CR(I,J,1-12)            
      DIMENSION QF(IB+2,JB+1,3),QC(IB+1,JB,3),DS(IB+1,JB,10)                     
      DIMENSION CR(IB+1,JB,12) 
                              !normale e tangenti a pannello generico co
      DIMENSION DSS(10) 
      COMMON/NO1/ CROOT,CTIP,XTIP,ZTIP,B,S,AR,PAY 
!                                                                       
!     B -IS FULL SPAN, C -ROOT CHORD, S - AREA                          
    S=0.5*B*(CROOT+CTIP) 
    C=S/B 
    AR=B*B/S 
!                                                                       
    IB1 = IB +1     
    JB1 = JB+1 
    IB2 = IB+2 
!                                                                       
!     CALCULATE PANEL CORNERPOINTS;  QF(I,J,(X,Y,Z))                    
!                                                                       
    CALL BODY(IB,JB,QF,thic)
      DO J = 1,JB1 
!mm        Y = -B/2.*cos(PAY*FLOAT(J-1)/JB) 
	y = -B/2.*cos(pay*float(j-1)/jb) 
        Yadim = 2.0*Y/B 
!       B - FULL SPAN, DXLE - LOCAL SWEEP                               
        DXLE = XTIP*Yadim 
        DZLE = ZTIP*Yadim 
        CHORD = CROOT-(CROOT-CTIP)*Yadim 
        DO I = 1,IB1 
          QF(I,J,1) = QF(I,1,1)*CHORD+DXLE 
          QF(I,J,2) = Y 
          QF(I,J,3) = QF(I,1,3)*CHORD+DZLE 
!         IF(J.EQ.1)write(23,*)I,QF(I,J,1),QF(I,J,2),QF(I,J,3),'I QF'   
              !I                                                        
        ENDDO 
!                                                                       
!       WAKE SHED PANEL POINTS (QF - IS IN BODY FRAME OF REFERENCE)     
        QF(IB2,J,1) = QF(IB1,J,1)+DXW 
        QF(IB2,J,2) = QF(IB1,J,2) 
                                  !+DZW                                 
        QF(IB2,J,3) = QF(IB1,J,3) 
            !J                                                          
      ENDDO
! -------------------------------------
!   ROTAZIONE DEL PROFILO DI ANGOLO ALFA e TRASLAZIONE HEAVE
! -------------------------------------
DO I = 1,IB2
    DO J = 1,JB1
	QF1 = QF(I,J,1)
	QF3 = QF(I,J,3)
	QF(I,J,1) =  QF1*CSA + QF3*SNA
	QF(I,J,3) = -QF1*SNA + QF3*CSA + SZ
    END DO
END DO     
!                                                                       
!     WING COLLOCATION POINTS                                           
!                                                                       
      DO J=1,JB 
        DO I=1,IB1 
          QC(I,J,1)=(QF(I,J,1)+QF(I,J+1,1)+QF(I+1,J+1,1)+QF(I+1,J,1))/4 
          QC(I,J,2)=(QF(I,J,2)+QF(I,J+1,2)+QF(I+1,J+1,2)+QF(I+1,J,2))/4 
          QC(I,J,3)=(QF(I,J,3)+QF(I,J+1,3)+QF(I+1,J+1,3)+QF(I+1,J,3))/4 
        ENDDO 
      ENDDO 
!  la geometria è stata lasciata nel SR corpo                          
!                                                                       
!     COMPUTATION OF TANGENTIAL CHORDWISE AND SPANWISE VECTORS          
!     DS(I,J,1-3) and DS(I,J,4-6)                                       
!     NORMAL VECTOR DS(I,J,7-9), PANEL AREA DS(I,J,10)                  
!     AND SOURCE STRENGTH (SIGMA)                                       
!                                                                       
      DO I=1,IB1 
        DO J=1,JB 
          CALL PANEL(                                                   &
     &    QF(I,J,1),QF(I,J,2),QF(I,J,3),                                &
     &    QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),                          &
     &    QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),                          &
     &    QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),                    &
     &    DS(I,J,1),DS(I,J,2),DS(I,J,3),DS(I,J,4),DS(I,J,5),DS(I,J,6),  &
     &    DS(I,J,7),DS(I,J,8),DS(I,J,9),DS(I,J,10))                     
        ENDDO 
      ENDDO 
!                                                                        
!                                                                       
!     TRANSFORM THE 4 PANEL CORNER POINTS INTO PANEL FRAME OF REF.      
!     THIS IS NEEDED LATER TO CALCULATE THE INFLUENCE COEFFICIENTS      
!                                                                       
!      i,j ------------ i,j+1     CR(i,j,1 - 3) = (X,Y,Z) corner I  ,J  
!          I1        4I           CR(i,j,4 - 6) = (X,Y,Z) corner I+1,J  
!          I          I           CR(i,j,7 - 9) = (X,Y,Z) corner I+1,J+1
!          I          I           CR(i,j,10-12) = (X,Y,Z) corner I  ,J+1
!          I2        3I                                                 
!    i+1,j ------------ i+1,j+1                                         
!                                                                       
!     TRANSFORMATION OF A FIELD POINT QC(I,J,1),QC(I,J,2),QC(I,J,3)     
!     INTO PANEL COORDINATES CR(I,J,1),CR(I,J,2),CR(I,J,3)              
!     QF(I,J,1),QF(I,J,2),QF(I,J,3) = PANEL CORNERPOINT                 
!     CR = coordinate nel SR ala con centro in QC                       
!                                                                       
      DO J=1,JB 
        DO I=1,IB1 
         DO K=1,10 
           DSS(K)=DS(I,J,K) 
         ENDDO 
!                                                                       
!        CALL CONVERT(QC(I,J,1),QC(I,J,2),QC(I,J,3),                    
!    &   QF(I,J,1),QF(I,J,2),QF(I,J,3),                                 
!    &   DSS,CR(I,J,1),CR(I,J,2),CR(I,J,3))                             
         X = QF(I,J,1) - QC(I,J,1) 
         Y = QF(I,J,2) - QC(I,J,2) 
         Z = QF(I,J,3) - QC(I,J,3) 
         CALL CONVERT1(X,Y,Z,DSS,CR(I,J,1),CR(I,J,2),CR(I,J,3)) 
!                                                                       
!        CALL CONVERT(QC(I,J,1),QC(I,J,2),QC(I,J,3),                    
!    &   QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),                           
!    &   DSS,CR(I,J,4),CR(I,J,5),CR(I,J,6))                             
         X = QF(I+1,J,1) - QC(I,J,1) 
         Y = QF(I+1,J,2) - QC(I,J,2) 
         Z = QF(I+1,J,3) - QC(I,J,3) 
         CALL CONVERT1(X,Y,Z,DSS,CR(I,J,4),CR(I,J,5),CR(I,J,6)) 
!                                                                       
!        CALL CONVERT(QC(I,J,1),QC(I,J,2),QC(I,J,3),                    
!    &   QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),                     
!    &   DSS,CR(I,J,7),CR(I,J,8),CR(I,J,9))                             
         X = QF(I+1,J+1,1) - QC(I,J,1) 
         Y = QF(I+1,J+1,2) - QC(I,J,2) 
         Z = QF(I+1,J+1,3) - QC(I,J,3) 
         CALL CONVERT1(X,Y,Z,DSS,CR(I,J,7),CR(I,J,8),CR(I,J,9)) 
!                                                                       
!        CALL CONVERT(QC(I,J,1),QC(I,J,2),QC(I,J,3),                    
!    &   QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),                           
!    &   DSS,CR(I,J,10),CR(I,J,11),CR(I,J,12))                          
         X = QF(I,J+1,1) - QC(I,J,1) 
         Y = QF(I,J+1,2) - QC(I,J,2) 
         Z = QF(I,J+1,3) - QC(I,J,3) 
         CALL CONVERT1(X,Y,Z,DSS,CR(I,J,10),CR(I,J,11),CR(I,J,12)) 
!                                                                       
              !I                                                        
        ENDDO 
            !J                                                          
      ENDDO 
!                                                                       
             !GRID                                                      
      RETURN 
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE PANEL(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,C1,C2,C3,    &
     & T1,T2,T3,V1,V2,V3,S)                                             
!     CALCOLA I VETTORI TANGENTE (x e y), NORMALE e AREA SU UN PANNELLO 
!     X,Y,Z-PANEL CORNERPOINTS, C,T,V-CHORDWISE, TANGENTIAL, NORMAL VECT
!     see K&P p. 251 e p 342                                            
!                                                                       
!          1--------3                                                   
!          I        I                                                   
!          I        I                                                   
!          I        I                                                   
!          2--------4                                                   
!                                                                       
!     FIRST CALCULATE CHORDWISE VECTOR tx                               
      A1=((X2+X4)-(X1+X3))/2.0 
      A2=((Y2+Y4)-(Y1+Y3))/2.0 
      A3=((Z2+Z4)-(Z1+Z3))/2.0 
      AA=SQRT(A1**2+A2**2+A3**2) 
      C1=A1/AA 
      C2=A2/AA 
      C3=A3/AA 
!                                                                       
!     NEXT, ANOTHER VECTOR IN THIS PLANE                                
      B1=X4-X1 
      B2=Y4-Y1 
      B3=Z4-Z1 
!                                                                       
!     NORMAL VECTOR                                                     
      V1=C2*B3-C3*B2 
      V2=B1*C3-C1*B3 
      V3=C1*B2-C2*B1 
      VV=SQRT(V1**2+V2**2+V3**2) 
      V1=V1/VV 
      V2=V2/VV 
      V3=V3/VV 
!                                                                       
!     TANGENTIAL VECTOR                                                 
      T1=V2*C3-V3*C2 
      T2=C1*V3-V1*C3 
      T3=V1*C2-V2*C1 
!                                                                       
!     CALCULATION OF PANEL AREA                                         
      E1=X3-X1 
      E2=Y3-Y1 
      E3=Z3-Z1 
      F1=X2-X1 
      F2=Y2-Y1 
      F3=Z2-Z1 
!                                                                       
!     NORMAL AREAS (F*B+B*E)                                            
      S11=F2*B3-F3*B2 
      S12=B1*F3-F1*B3 
      S13=F1*B2-F2*B1 
      S21=B2*E3-B3*E2 
      S22=E1*B3-B1*E3 
      S23=B1*E2-B2*E1 
      S=0.5*(SQRT(S11**2+S12**2+S13**2)+SQRT(S21**2+S22**2+S23**2)) 
      RETURN 
          !PANEL                                                        
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE WING(IB,JB,QC,CR,I,J,DDUBJ,DS,A1,SIGMA,RH) 
!                                                                       
!     CALCULATES WING POTENTIAL AT A COLLOCATION POINT (I,J)            
!     DUE TO WING DOUBLET DISTRIBUTION ON PANEL (I1,J1)                 
!     IN A WING FIXED COORDINATE SYSTEM                                 
!     I,J = COLLOCATION POINT, I1,J1 = INFLUENCING PANEL                
!     see K&P par 12.5 pag 353 sgg                                      
!                                                                       
                              !coordinate centro pannelli corpo         
      DIMENSION QC(IB+1,JB,3) 
                               !coordinate 4 vertici pannelli corpo     
      DIMENSION CR(IB+1,JB,12) 
      DIMENSION DDUBJ(JB),A1(IB,JB) 
                               !normale e tangenti a pannello corpo     
      DIMENSION DS(IB+1,JB,10) 
                             !sorgenti su corpo                         
      DIMENSION SIGMA(IB,JB) 
                          ! vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                              !normale e tangenti a pannello generico co
      DIMENSION DSS(10) 
!                                                                       
!     WING FIXED PANEL CORNERPOINTS                                     
!     CR(I,J,1)=X(i,j); CR(I,J,2)=Y(I,J); CR(I,J,3)=Z(I,J)              
! calcola il potenziale in P(x,y,z) prodotta dal pannello I,J           
!     i,j I---1----I i,j+1                                              
!         I        I                                                    
!         4  I,J   2                         P(x,y,z)                   
!         I        I                                                    
!   i+1,j I---3----I i+1,j+1                                            
!                                                                       
      IB1=IB+1 
      RH=0. 
!     INFLUENCING PANEL I1,J1                                           
      DO I1=1,IB 
        DO J1=1,JB 
!         CALCULATE WAKE CONTRIBUTION (effetto pannelli scia vicina i=IB
                                    !pannello al TE ventre              
          IF(I1.EQ.1) THEN 
!           FIRST CONVERT COLLOCATION POINT (IB1,J1) TO PANEL COORDINATE
!           AND THEN CALCULATE INFLUENCE COEFFICIENTS                   
            DO K=1,10 
              DSS(K)=DS(IB1,J1,K) 
            ENDDO 
!           CALL CONVERT(QC(IB1,J1,1),QC(IB1,J1,2),QC(IB1,J1,3),        
!    &                   QC(I,J,1),QC(I,J,2),QC(I,J,3),                 
!    &                   DSS,XC,YC,ZC)                                  
            X = QC(I,J,1) - QC(IB1,J1,1) 
            Y = QC(I,J,2) - QC(IB1,J1,2) 
            Z = QC(I,J,3) - QC(IB1,J1,3) 
            CALL CONVERT1(X,Y,Z,DSS,XC,YC,ZC) 
!                                                                       
            DO K=1,12 
              CRR(K)=CR(IB1,J1,K) 
            ENDDO 
            CALL INFLUENCE(DMU2,DSIG,XC,YC,ZC,CRR) 
!           write(*,*)i,j,j1,DMU2,' i,j,j1,DMU2' !cfr con calc semplif  
!                                                                       
! riv segno giusto?                                                     
            DDUBJ(J1)=-DMU2 
          ELSE 
            DMU2=0.0 
          ENDIF 
                                       !pannello al TE dorso            
          IF(I1.EQ.IB) DMU2=+DDUBJ(J1) 
!         END OF WAKE INFLUENCE CALCULATION                             
!                                                                       
          DO K=1,10 
            DSS(K)=DS(I1,J1,K) 
          ENDDO 
!         CONVERT COLLOCATION POINT QC(I,J,k) TO PANEL (I1,J1) COORDINAT
!  primi 3= coordin centro influencing panel (I1,J1),                   
!  secondi 3= colloc pt (I,J) da trasformare;                           
!  XC,YC,ZC= colloc pt nel sr pannello                                  
!         CALL CONVERT(QC(I1,J1,1),QC(I1,J1,2),QC(I1,J1,3),             
!    &                 QC(I,J,1),QC(I,J,2),QC(I,J,3),                   
!    &                 DSS,XC,YC,ZC)                                    
          X = QC(I,J,1) - QC(I1,J1,1) 
          Y = QC(I,J,2) - QC(I1,J1,2) 
          Z = QC(I,J,3) - QC(I1,J1,3) 
          CALL CONVERT1(X,Y,Z,DSS,XC,YC,ZC) 
!                                                                       
!         effetto del pannello I1,J1 sul collocation pt XC              
          DO K=1,12 
            CRR(K)=CR(I1,J1,K) 
          ENDDO 
          CALL INFLUENCE(DMU,DSIG,XC,YC,ZC,CRR) 
!         CALL P12(DMU,DSIG,XC,YC,ZC,CRR)                               
!         CALL INFLUENCEFF(DMUFF,DSIGFF,XC,YC,ZC,CRR,                   
!    &                       VXS,VYS,VZS,VXD,VYD,VZD,dd)                
! calcolo del ctb FAR-FIELD della sorgente nel SR ground                
!         R0=SQRT((QC(I,J,1)-QC(I1,J1,1))**2+                           
!    &        (QC(I,J,2)-QC(I1,J1,2))**2+(QC(I,J,3)-QC(I1,J1,3))**2)    
!         SOU_FF=DSS(10)/(4.*R0*ACOS(-1.))                              
!         write(128,10)I,J,I1,J1,DMUFF,DMU,DSIGFF,DSIG,dd               
!  10 FORMAT(4(I2,1x),5(F16.8,1x))                                      
!                                                                       
!         A PANEL INFLUENCE ON ITSELF IS DMU=1/2                        
          IF((I1.EQ.I).AND.(J1.EQ.J)) DMU=+0.5 
!                                                                       
!         A1(I1,J1) - IS THE INFLUENCE MATRIX COEFFICIENT               
          A1(I1,J1) = DMU -DMU2 
          RH = RH-DSIG*SIGMA(I1,J1) 
!                                                                       
               !J1                                                      
        ENDDO 
               !I1                                                      
      ENDDO 
!                                                                       
      RETURN 
          !WING                                                         
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE CONVERT1(XB,YB,ZB,DSS,XP,YP,ZP) 
!     DIMENSION DS(IB1,JB,10)                                           
      DIMENSION DSS(10) 
!                                                                       
!     TRANSFORMATION OF A FIELD POINT XO,YO,ZO                          
!     INTO PANEL COORDINATES XP,YP,ZP                                   
!     XO = coordin centro influencing panel,                            
!     XB colloc pt da trasformare (PANEL  CORNERPOINT)                  
!     XP= colloc pt nel sr pannello                                     
!     DSS(1-3),DSS(4-6),DSS(/-9) ARE CHORDWISE, TANGENTIAL, AND NORMAL V
!                                                                       
      XP=XB*DSS(1)+YB*DSS(2)+ZB*DSS(3) 
      YP=XB*DSS(4)+YB*DSS(5)+ZB*DSS(6) 
      ZP=XB*DSS(7)+YB*DSS(8)+ZB*DSS(9) 
!                                                                       
      RETURN 
           !CONVERT1                                                    
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE INFLUENCE(DOUB,SOUR,XC,YC,ZC,CRR) 
      DIMENSION CRR(12) 
      COMMON/NO3/ FOURPAY 
!                                                                       
!     DOUBLET AND SOURCE POTENTIAL AT POINT (XC,YC,ZC) DUE TO           
!     PANEL I,J  i.e. CR(I,J,1-2-3) (X1,Y1,Z1,...X4,Y4,Z4),             
!     SEE KATZ & PLOTKIN PP 283-6                                       
!     SEE ALSO Hess & Smith 1967 p 49-58                                
!     SEE ALSO Ortega et al. 2010 p 11                                  
!     SEE ALSO Maskew 1987 p 37                                         
!                                                                       
      ONE=1. 
      EPS=EPSILON(ONE)*10. 
      COE=1./FOURPAY 
!                                                                       
      X1 = CRR(1) 
      Y1 = CRR(2) 
      Z1 = CRR(3) 
      X2 = CRR(4) 
      Y2 = CRR(5) 
      Z2 = CRR(6) 
      X3 = CRR(7) 
      Y3 = CRR(8) 
      Z3 = CRR(9) 
      X4 = CRR(10) 
      Y4 = CRR(11) 
      Z4 = CRR(12) 
!                                                                       
      if(abs(z1)+abs(z2)+abs(z3)+abs(z4) .GE. EPS) THEN 
        write(*,*)z1,z2,z3,z4,'z 632' 
        write(*,*)xc,xy,zc,'xc,xy,zc' 
        stop '768 sub INFLUENCE' 
      ENDIF 
!                                                                       
      PNLX=.25*(X1+X2+X3+X4)                                            
      PNLY=.25*(Y1+Y2+Y3+Y4)                                            
      PNLZ=.25*(Z1+Z2+Z3+Z4)                                            
      PNX=XC-PNLX                                                       
      PNY=YC-PNLY                                                       
      PNZ=ZC-PNLZ                                                       
!     PNLX è X0 ecc PNS2 è R0q, PNS è R0                             
      PNS2=PNX*PNX+PNY*PNY+PNZ*PNZ                                      
      PNS=SQRT(PNS2)                                                    
      PNS3=PNS*PNS2                                                     
!                                                                       
!                                                                       
! provvis x cfr calcolo area pannello                                   
      D1X=X3-X1                                                         
!     D1Y=Y3-Y1                                                         
!     D1Z=Z3-Z1                                                         
!     D2X=X4-X2                                                         
      D2Y=Y4-Y2                                                         
!     D2Z=Z4-Z2                                                         
!     CRX=D1Y*D2Z-D2Y*D1Z                                               
!     CRY=D2X*D1Z-D1X*D2Z                                               
!     CRZ=D1X*D2Y-D2X*D1Y                                               
!     CRSQ=SQRT(CRX*CRX+CRY*CRY+CRZ*CRZ)                                
!     AREA1=CRSQ/2.                                                     
!                                                                       
      AREA=D1X*D2Y                                                      
!     DD=PNS/SQRT(ABS(AREA))                                            
!     RNX=CRX/CRSQ                                                      
!     RNY=CRY/CRSQ                                                      
!     RNZ=CRZ/CRSQ                                                      
!     write(*,*)PNX*RNX+PNY*RNY+PNZ*RNZ,PNZ !PN dot N                   
      R23x=0.5*(X2+X3)-PNLX                                             
      R23y=0.5*(Y2+Y3)-PNLY                                             
      R23z=0.5*(Z2+Z3)-PNLZ                                             
      SMP=SQRT(R23x**2+R23y**2+R23z**2)                                 
      R34x=0.5*(X3+X4)-PNLX                                             
      R34y=0.5*(Y3+Y4)-PNLY                                             
      R34z=0.5*(Z3+Z4)-PNLZ                                             
      SMQ=SQRT(R34x**2+R34y**2+R34z**2)                                 
      AVGSIZE=SMP+SMQ                                                   
      DD=AVGSIZE                                                        
!     write(*,*)AVGSIZE,' DIMENSIONE PANNELLO',PNS,' DISTANZA PTO'      
!                                                                       
!     PNLX è X0 ecc PNS2 è R0q, PNS è R0                             
!     T1=AREA*COE                                                       
!     DD=PNS/SQRT(ABS(AREA))                                            
                                                                        
      DOUB= 0.                                                          
      SOUR= 0.0                                                         
!                                                                       
!     SOURCE (S, SOUR) AND DOUBLET (Q, DOUB) INFLUENCE IN PANEL COORDINA
!     FOR TRIANGULAR PANEL THE 4TH SIDE CONTRIBUTION IS ZERO            
!                                                                       
!     calcolo FAR-FIELD x doppietta e sorgente                          
!     FF ok x doppietta se DD > 10 mentre x sorg > 40                   
!     DOUBLET CONTRIBUTION K&P p 249                                    
      IF(DD.GT.50.AND.ABS(PNS3).GT.EPS) THEN                            
        write(*,*)'FAR-FIELD SORG+DOPPIETTA'                            
        T1=AREA*COE                                                     
        DOUB=-T1*PNZ/PNS3                                               
        SOUR= T1/PNS                                                    
        RETURN                                                          
      ENDIF                                                             
      IFF=0                                                             
      IF(DD.GT.20.AND.ABS(PNS3).GT.EPS) IFF=1                           
!                                                                       
!     EPS, PANEL SIDE CUTOFF DISTANCE                                   
!     PANEL SIDE (D) DISTANCE (R), E, AND H (EQS. 10.90 & 10.92-10.94)  
!                                                                       
      X2mX1=X2-X1                                                       
      X3mX2=X3-X2                                                       
      X4mX3=X4-X3                                                       
      X1mX4=X1-X4                                                       
      Y2mY1=Y2-Y1                                                       
      Y3mY2=Y3-Y2                                                       
      Y4mY3=Y4-Y3                                                       
      Y1mY4=Y1-Y4                                                       
      XCmX1=XC-X1                                                       
      XCmX2=XC-X2                                                       
      XCmX3=XC-X3                                                       
      XCmX4=XC-X4                                                       
      YCmY1=YC-Y1                                                       
      YCmY2=YC-Y2                                                       
      YCmY3=YC-Y3                                                       
      YCmY4=YC-Y4                                                       
!                                                                       
      X2mX1q=X2mX1**2                                                   
      X3mX2q=X3mX2**2                                                   
      X4mX3q=X4mX3**2                                                   
      X1mX4q=X1mX4**2                                                   
      Y2mY1q=Y2mY1**2                                                   
      Y3mY2q=Y3mY2**2                                                   
      Y4mY3q=Y4mY3**2                                                   
      Y1mY4q=Y1mY4**2                                                   
      XCmX1q=XCmX1**2                                                   
      XCmX2q=XCmX2**2                                                   
      XCmX3q=XCmX3**2                                                   
      XCmX4q=XCmX4**2                                                   
      YCmY1q=YCmY1**2                                                   
      YCmY2q=YCmY2**2                                                   
      YCmY3q=YCmY3**2                                                   
      YCmY4q=YCmY4**2                                                   
      ZCq=ZC*ZC                                                         
!                                                                       
      D1=SQRT(X2mX1q+Y2mY1q)                                            
      D2=SQRT(X3mX2q+Y3mY2q)                                            
      D3=SQRT(X4mX3q+Y4mY3q)                                            
      D4=SQRT(X1mX4q+Y1mY4q)                                            
!                                                                       
      R1=SQRT(XCmX1q+YCmY1q+ZCq)                                        
      R2=SQRT(XCmX2q+YCmY2q+ZCq)                                        
      R3=SQRT(XCmX3q+YCmY3q+ZCq)                                        
      R4=SQRT(XCmX4q+YCmY4q+ZCq)                                        
!                                                                       
!     SOURCE CONTRIBUTION K&P p 247 H&S p 57                            
!                                                                       
      S1=0.                                                             
      TT=R1+R2+D1                                                       
      IF(TT.GT.EPS.AND.D1.GT.EPS)THEN                                   
!       R2=B; R1=A; ZC=PN; X2-X1=SM; D1=S                               
        S1=-(XCmX1*Y2mY1-YCmY1*X2mX1)/D1*LOG(TT/(R1+R2-D1))                                          
      ENDIF                                                             
!                                                                       
      S2=0.                                                             
      TT=R2+R3+D2                                                       
      IF(TT.GT.EPS.AND.D2.GT.EPS)THEN                                   
        S2=-(XCmX2*Y3mY2-YCmY2*X3mX2)/D2*LOG(TT/(R2+R3-D2))                                          
      ENDIF                                                             
!                                                                       
      S3=0.                                                             
      TT=R3+R4+D3                                                       
      IF(TT.GT.EPS.AND.D3.GT.EPS)THEN                                   
        S3=-(XCmX3*Y4mY3-YCmY3*X4mX3)/D3*LOG(TT/(R3+R4-D3))                                         
      ENDIF                                                             
!                                                                       
      S4=0.                                                             
      TT=R4+R1+D4                                                       
      IF(TT.GT.EPS.AND.D4.GT.EPS)THEN                                   
        S4=-(XCmX4*Y1mY4-YCmY4*X1mX4)/D4*LOG(TT/(R4+R1-D4))                                         
      ENDIF                                                             
!                                                                       
!     ADD CONTRIBUTIONS FROM THE 4 SIDES                                
!                                                                       
      SOUR1= (S1+S2+S3+S4)*COE                                          
!                                                                       
      IF(IFF.EQ.1) THEN                                                 
!       calcolo FAR-FIELD x DOPPIETTA                                   
        write(*,*)'--- FAR-FIELD DOPPIETTA'                             
        DOUB=-AREA*PNZ/PNS3                                             
!                   NON con ABS(ZC)                                     
        SOUR= SOUR1 +    (ZC)*DOUB                                      
        RETURN                                                          
      ENDIF                                                             
!                                                                       
!     DOUBLET CONTRIBUTION                                              
                                                                        
      IF(ABS(ZC).LT.EPS) GOTO 100                                       
!                                                                       
      E1= XCmX1q+ZCq                                                    
      E2= XCmX2q+ZCq                                                    
      E3= XCmX3q+ZCq                                                    
      E4= XCmX4q+ZCq                                                    
!                                                                       
      H1= XCmX1 * YCmY1                                                 
      H2= XCmX2 * YCmY2                                                 
      H3= XCmX3 * YCmY3                                                 
      H4= XCmX4 * YCmY4                                                 
!                                                                       
      Q1=0.                                                             
      IF (D1.GT.EPS) THEN                                               
        F=Y2mY1*E1-X2mX1*H1                                             
        G=Y2mY1*E2-X2mX1*H2                                             
!       R2=B; R1=A; ZC=PN; X2-X1=SM; D1=S                               
        RNUM=ZC*X2mX1*(F*R2-G*R1)                                       
        DNOM=ZCq*X2mX1q*R1*R2 + F*G                                     
        Q1=ATAN2(RNUM,DNOM)                                             
      ENDIF                                                             
!                                                                       
      Q2=0.                                                             
      IF (D2.GT.EPS) THEN                                               
        F=Y3mY2*E2-X3mX2*H2                                             
        G=Y3mY2*E3-X3mX2*H3                                             
        RNUM=ZC*X3mX2*(F*R3-G*R2)                                       
        DNOM=ZCq*X3mX2q*R2*R3 + F*G                                     
        Q2=ATAN2(RNUM,DNOM)                                             
      ENDIF                                                             
!                                                                       
      Q3=0.                                                             
      IF (D3.GT.EPS) THEN                                               
        F=Y4mY3*E3-X4mX3*H3                                             
        G=Y4mY3*E4-X4mX3*H4                                             
        RNUM=ZC*X4mX3*(F*R4-G*R3)                                       
        DNOM=ZCq*X4mX3q*R3*R4 + F*G                                     
        Q3=ATAN2(RNUM,DNOM)                                             
      ENDIF                                                             
!                                                                       
      Q4=0.                                                             
      IF (D4.GT.EPS) THEN                                               
        F=Y1mY4*E4-X1mX4*H4                                             
        G=Y1mY4*E1-X1mX4*H1                                             
        RNUM=ZC*X1mX4*(F*R1-G*R4)                                       
        DNOM=ZCq*X1mX4q*R4*R1 + F*G                                     
        Q4=ATAN2(RNUM,DNOM)                                             
      ENDIF                                                             
!                                                                       
!     ADD CONTRIBUTIONS FROM THE 4 SIDES                                
      DOUB= (Q1+Q2+Q3+Q4)*COE                                           
!                                                                       
  100 CONTINUE                                                          
!                                                                       
!     SOUR= SOUR1 + ABS(ZC)*DOUB  !times source strength   10.89 NON con
      SOUR= SOUR1 + ZC*DOUB                                             
!                                                                       
      RETURN                                                            
                                        !semi-diagonal                  
                                        !semi-diagonal                  
                                              !FAR-FIELD solo x DOPPIETT
                              !d12 (10.90a)                             
                              !d23 (10.90a)                             
                              !d34 (10.90c)                             
                              !d41 (10.90d)                             
                                   !(10.92)                             
                                   !(10.92)                             
                                   !(10.92)                             
                                   !(10.92)                             
                                                    !(parte1 in 10.89)  
                                                    !(parte2 in 10.89)  
                                                    !(parte3 in 10.89)  
                                                    !(parte4 in 10.89)  
                                    !ctb parziale source 10.89          
                                    !times source strength   10.89      
                      !(10.93)                                          
                      !(10.93)                                          
                      !(10.93)                                          
                      !(10.93)                                          
                           !(10.94)                                     
                           !(10.94)                                     
                           !(10.94)                                     
                           !(10.94)                                     
                            !PA                                         
                            !PB                                         
                            !(10.105)                                   
                            !(10.105)                                   
                            !(10.105)                                   
                            !(10.105)                                   
                                             !times doublet strength (10
                             !times source strength   10.89             
          !INFLUENCE                                                    
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE DECOMP(N,NDIM,A,IP) 
      REAL A(NDIM,NDIM),T 
      INTEGER IP(NDIM) 
!     MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.                 
!     N = ORDER OF MATRIX. NDIM = DECLARED DIMENSION OF ARRAY A.        
!     A = MATRIX TO BE TRIANGULARIZED.                                  
!     IP(K) , K .LT. N = INDEX OF K-TH PIVOT ROW.                       
!                                                                       
!     EPS = 0.000001                                                    
      ONE=1. 
      EPS=EPSILON(ONE)*10. 
!                                                                       
      IP(N) = 1 
      DO K = 1, N 
        IF(K.EQ.N) GOTO 5 
        KP1 = K + 1 
        M=K 
        DO I = KP1, N 
          IF( ABS(A(I,K)).GT.ABS(A(M,K))) M=I 
        ENDDO 
        IP(K) = M 
        IF(M.NE.K) IP(N) = -IP(N) 
        T = A(M,K) 
        A(M,K) = A(K,K) 
        A(K,K) = T 
        IF(abs(T) .LE. EPS) GO TO 5 
        DO I = KP1, N 
          A(I,K) = -A(I,K)/T 
        ENDDO 
        DO J = KP1, N 
          T = A(M,J) 
          A(M,J) = A(K,J) 
          A(K,J) = T 
          IF(abs(T) .LE. EPS) GO TO 4 
          DO I = KP1, N 
            A(I,J) = A(I,J) + A(I,K)*T 
          ENDDO 
    4   CONTINUE 
        ENDDO 
    5   IF(abs(A(K,K)) .LE. EPS) IP(N) = 0 
      ENDDO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE SOLVER(N,NDIM,A,B,IP) 
      REAL A(NDIM,NDIM), B(NDIM), T 
      INTEGER IP(NDIM) 
!     SOLUTION OF LINEAR SYSTEM, A*X = B.                               
!     N = ORDER OF MATRIX.                                              
!     NDIM = DECLARED DIMENSION OF THE ARRAY A.                         
!     B = RIGHT HAND SIDE VECTOR.                                       
!     IP = PIVOT VECTOR OBTAINED FROM SUBROUTINE DECOMP.                
!     B = SOLUTION VECTOR, X.                                           
!                                                                       
      IF(N.EQ.1) GOTO 9 
      NM1 = N - 1 
      DO K = 1, NM1 
        KP1 = K + 1 
        M = IP(K) 
        T = B(M) 
        B(M) = B(K) 
        B(K) = T 
        DO I = KP1, N 
          B(I) = B(I) + A(I,K)*T 
        ENDDO 
      ENDDO 
      DO KB = 1, NM1 
        KM1 = N - KB 
        K = KM1 + 1 
        B(K) = B(K)/A(K,K) 
        T = -B(K) 
        DO I = 1, KM1 
          B(I) = B(I) + A(I,K)*T 
        ENDDO 
      ENDDO 
    9 B(1) = B(1)/A(1,1) 
      RETURN 
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE POT_WAKE(IB,JB,NSTEPS,QFW,QC,DUBW,I,J,N,RH) 
!                                                                       
!     CALCULATES POTENTIAL AT A COLLOCATION POINT (I,J) on WING         
!     DUE TO WAKE DOUBLET DISTRIBUTION DUBW(N,JJ) AT T=N*DT             
!     IN A WING FIXED COORDINATE SYSTEM                                 
!                                                                       
                          ! centri pannelli scia                        
      DIMENSION QCW(3) 
                                      !corner points inizio pannelli sci
      DIMENSION QFW(NSTEPS+1,JB+1,3) 
                              !coordinate centro pannelli corpo         
      DIMENSION QC(IB+1,JB,3) 
                                !intensità doppiette scia              
      DIMENSION DUBW(NSTEPS,JB) 
                          ! vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                              !normale e tangenti a pannello generico co
      DIMENSION DSS(10) 
!                                                                       
      IB1=IB+1 
      RH=0. 
      DO IT=1,N-1 
        DO JJ=1,JB 
!         WAKE PANEL CENTER (origine of the local R.S.)                 
          DO K=1,3 
            QCW(K)=(QFW(IT+1,JJ,K)+QFW(IT+1,JJ+1,K)+                    &
     &              QFW(IT,JJ,K)+QFW(IT,JJ+1,K))/4.                     
          ENDDO 
!                                                                       
!     COMPUTATION OF TANGENTIAL CHORDWISE AND SPANWISE VECTORS          
!     DSS(1-3) and DSS(4-6) NORMAL VECTOR DSS(7-9), PANEL AREA DSS(10)  
!     calcola il DSS del pannello IT,JJ                                 
!                                                                       
          CALL PANEL(                                                   &
     &    QFW(IT+1,JJ,1),QFW(IT+1,JJ,2),QFW(IT+1,JJ,3),                 &
     &    QFW(IT,JJ,1),QFW(IT,JJ,2),QFW(IT,JJ,3),                       &
     &    QFW(IT+1,JJ+1,1),QFW(IT+1,JJ+1,2),QFW(IT+1,JJ+1,3),           &
     &    QFW(IT,JJ+1,1),QFW(IT,JJ+1,2),QFW(IT,JJ+1,3),                 &
     &    DSS(1),DSS(2),DSS(3),DSS(4),DSS(5),                           &
     &    DSS(6),DSS(7),DSS(8),DSS(9),DSS(10))                          
!                                                                       
!         CONVERT COLLOCATION POINT TO WAKE PANEL COORDINATES           
!         il S.R. ha origine al centro del pannello di scia             
          XB=QC(I,J,1)-QCW(1) 
          YB=QC(I,J,2)-QCW(2) 
          ZB=QC(I,J,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,XC,YC,ZC) 
!                                                                       
!         TRANSFORMATION OF A WAKE POINT QFW(..,1/2/3)                  
!         INTO PANEL COORDINATES CRR(1),CRR(2),CRR(3)                   
!         QFW(I,J,1),QFW(I,J,2),QFW(I,J,3) = WAKE PANEL CORNERPOINT     
!         CRR = coordinate nel SR ala con centro in QCW                 
!                                                                       
          XB=QFW(IT+1,JJ,1)-QCW(1) 
          YB=QFW(IT+1,JJ,2)-QCW(2) 
          ZB=QFW(IT+1,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(1),CRR(2),CRR(3)) 
!                                                                       
          XB=QFW(IT,JJ,1)-QCW(1) 
          YB=QFW(IT,JJ,2)-QCW(2) 
          ZB=QFW(IT,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(4),CRR(5),CRR(6)) 
!                                                                       
          XB=QFW(IT,JJ+1,1)-QCW(1) 
          YB=QFW(IT,JJ+1,2)-QCW(2) 
          ZB=QFW(IT,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(7),CRR(8),CRR(9)) 
!                                                                       
          XB=QFW(IT+1,JJ+1,1)-QCW(1) 
          YB=QFW(IT+1,JJ+1,2)-QCW(2) 
          ZB=QFW(IT+1,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(10),CRR(11),CRR(12)) 
!                                                                       
          CALL INFLUENCEW(DMU,XC,YC,ZC,CRR) 
!         CALL INFLUENCEWFF(DMU,XC,YC,ZC,CRR)                           
!         CALL P12W(DMU,XC,YC,ZC,CRR)                                   
!                                                                       
          RH = RH + DMU*DUBW(IT,JJ) 
!                                                                       
               !JJ                                                      
        ENDDO 
               !IT                                                      
      ENDDO 
!                                                                       
      RETURN 
          !POT_WAKE                                                     
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE FORCES(IB,JB,QF,QC,DS,DUB,DUBPREC,FS,FSPREC,DT,UBC,WBC,N,DL,DD,CP,CL,CD,CM,FL,QINF,QL,QM) 
      DIMENSION QF(IB+2,JB+1,3),QC(IB+1,JB,3),DS(IB+1,JB,10),UBC(:),WBC(:)
      DIMENSION DUB(IB,JB),dubprec(ib,jb),FS(IB,JB),FSPREC(IB,JB),DL(IB,JB),DD(IB,JB),CP(IB,JB)
	integer:: n 
      COMMON/NO1/ CROOT,CTIP,XTIP,ZTIP,B,S,AR,PAY 
      COMMON/NO2/ RO,VT,UT 
!                                                                       
!     ==================                                                
!     FORCES CALCULATION                                                
!     ==================                                                
!                                                                       
      FL=0. 
      FD=0. 
      FM=0. 
      QUE=0.5*RO*UT*UT 
    IM = (IB+2)/2
!                                                                       
    open(14,file='qinflmcp.txt')
    write(14,*) 'I  J	QI  QM	QL  CP	DS7	DS9'
      DO J=1,JB 
        DO I=1,IB 
          I1=I-1 
          I2=I+1 
          J1=J-1 
          J2=J+1 
          IF(I.EQ.1) I1=1 
          IF(I.EQ.IB) I2=IB 
          IF(J.EQ.1) J1=1 
          IF(J.EQ.JB) J2=JB 
!         CHORDWISE VELOCITY                                            
          XF=0.5*(QF(I+1,J,1)+QF(I+1,J+1,1)) 
          YF=0.5*(QF(I+1,J,2)+QF(I+1,J+1,2)) 
          ZF=0.5*(QF(I+1,J,3)+QF(I+1,J+1,3)) 
          XR=0.5*(QF(I,J,1)+QF(I,J+1,1)) 
          YR=0.5*(QF(I,J,2)+QF(I,J+1,2)) 
          ZR=0.5*(QF(I,J,3)+QF(I,J+1,3)) 
          DX2=QC(I2,J,1)-XF 
          DY2=QC(I2,J,2)-YF 
          DZ2=QC(I2,J,3)-ZF 
          DX3=QC(I1,J,1)-XR 
          DY3=QC(I1,J,2)-YR 
          DZ3=QC(I1,J,3)-ZR 
          DL1=SQRT((XF-XR)**2+(YF-YR)**2+(ZF-ZR)**2) 
          DL2=SQRT(DX2**2+DY2**2+DZ2**2) 
          DL3=SQRT(DX3**2+DY3**2+DZ3**2) 
          DLL=DL1+DL2+DL3 
          IF(I.EQ.1) DLL=DL1/2.0+DL2 
          IF(I.EQ.IB) DLL=DL1/2.0+DL3 
          QL=-(DUB(I2,J)-DUB(I1,J))/DLL 
!                                                                       
!         SPANWISE VELOCITY                                             
!                                                                       
          DX=QC(I,J2,1)-QC(I,J1,1) 
          DY=QC(I,J2,2)-QC(I,J1,2) 
          DZ=QC(I,J2,3)-QC(I,J1,3) 
          DR=SQRT(DX**2+DY**2+DZ**2) 
          QM=-(DUB(I,J2)-DUB(I,J1))/DR 
!                                                                       
!         FIRST ORDER CORRECTION FOR PANEL SWEEP                        
!                                                                       
          QL=QL+QM*(DX**2+DZ**2)/DR 
          QM=QM*(DY**2+DZ**2)/DR 
          QINF= -UBC(I)*DS(I,J,1) - WBC(I)*DS(I,J,3) 
          !capire
if (n.eq.1) dubprec = 0.

	CP(I,J)=1.0-((QINF+QL)**2+QM**2)/(UT**2) - 2*(FS(I,J)-FSPREC(I,J))/(DT*UT**2)
                                                  !=-CP*Area*Nz         
          DL(I,J)=-CP(I,J)*DS(I,J,10)*DS(I,J,9) 
                                                  != CP*Area*nx         
          DD(I,J)=CP(I,J)*DS(I,J,10)*DS(I,J,7) 
          FL=FL+DL(I,J) 
          FD=FD+DD(I,J) 
          FM=FM+DL(I,J)*QC(I,J,1) 
    write(14,*) i,j,qinf,qm,ql,cp(i,j),DS(i,j,7),DS(i,j,9)
        ENDDO 
      ENDDO 
      CL=FL/(QUE*S)/2. !tt correzione cl. serve anche a cm e cd???? 
      CD=FD/(QUE*S) 
      CM=FM/(QUE*S*CROOT) 
      RETURN 
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VELOCE2(IB,JB,QC,CR,DUB,SIGMA,DS,XX,YY,ZZ,U,V,W) 
!                                                                       
!     CALCULATES WAKE VELOCITY AT A COLLOCATION POINT XX,YY,ZZ i.e. QFW(
!     DUE TO WING SOURCE AND DOUBLET DISTRIBUTION ON PANEL (I1,J1)      
!     IN A WING FIXED COORDINATE SYSTEM                                 
!                                                                       
                               !coordinate centro pannelli corpo        
      DIMENSION QC(IB+1,JB,3) 
                               !coordinate 4 vertici pannelli corpo     
      DIMENSION CR(IB+1,JB,12) 
                             !doppiette su corpo                        
      DIMENSION DUB(IB,JB) 
                             !sorgenti su corpo                         
      DIMENSION SIGMA(IB,JB) 
                               !normale e tangenti a pannello corpo     
      DIMENSION DS(IB+1,JB,10) 
                           !vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                           !normale e tangenti a pannello generico corpo
      DIMENSION DSS(10) 
!                                                                       
      COMMON/NO3/ FOURPAY 
!                                                                       
!     WING FIXED PANEL CORNERPOINTS                                     
!     CR(I,J,1)=X(i,j); CR(I,J,2)=Y(I,J); CR(I,J,3)=Z(I,J)              
! calcola il potenziale in P(x,y,z) prodotta dal pannello I,J           
!     i,j I---1----I i,j+1                                              
!         I        I                                                    
!         4  I,J   2                         P(x,y,z)                   
!         I        I                                                    
!   i+1,j I---3----I i+1,j+1                                            
!                                                                       
      IB1=IB+1 
!                                                                       
!-riv X1=(X-X0)*CS1+(Z-Z0)*SN1                                          
!-riv Y1=Y                                                              
!-riv Z1=-(X-X0)*SN1+(Z-Z0)*CS1                                         
!     X,Y,Z è il punto della scia ove voglio calcolare la velocità    
!                                                                       
      U=0. 
      V=0. 
      W=0. 
!                                                                       
!     INFLUENCING PANEL I1,J1 sul corpo                                 
      DO I1=1,IB 
        DO J1=1,JB 
!                                                                       
          DO K=1,10 
            DSS(K)=DS(I1,J1,K) 
          ENDDO 
!         CONVERT COLLOCATION POINT XX,YY,ZZ TO PANEL (I1,J1) COORDINATE
!  primi 3= coordin centro influencing panel (I1,J1),                   
!  secondi 3= colloc pt (I,J) da trasformare; X,Y,Z= nodi della scia QFW
!  XC,YC,ZC= colloc pt nel sr pannello                                  
          X = XX - QC(I1,J1,1) 
          Y = YY - QC(I1,J1,2) 
          Z = ZZ - QC(I1,J1,3) 
                                                                        
          CALL CONVERT1(X,Y,Z,DSS,XC,YC,ZC) 
                                                                        
!         effetto del pannello I1,J1 sul collocation pt X,Y,Z           
          DO K=1,12 
            CRR(K)=CR(I1,J1,K) 
          ENDDO 
!                                                                       
!         CALL VEL_WING(CRR,XC,YC,ZC,DUB(I1,J1),SIGMA(I1,J1),           
!    &                  Utx,Uty,Un,UtxS,UtyS,UnS,UtxD,UtyD,UnD)         
! con valutazione FAR FIELD per pannelli lontani                        
          CALL VEL_WINGFF(CRR,XC,YC,ZC,DUB(I1,J1),SIGMA(I1,J1),         &
     &                  Utx,Uty,Un,UtxS,UtyS,UnS,UtxD,UtyD,UnD,PNS2)    
!                                                                       
!         proietta per calcolare le componenti nel S.R corpo            
!         il SR corpo ha origine nel LE e corda su asse x               
!         ux=utx*tx_x+uty*ty_x+un*n_x                                   
!         uy=utx*tx_y+uty*ty_y+un*n_y                                   
!         uz=utx*tx_z+uty*ty_z+un*n_z                                   
!                                                                       
!         VXS=UtxS*DS(I1,J1,1)+UtyS*DS(I1,J1,4)+UnS*DS(I1,J1,7)         
!         VYS=UtxS*DS(I1,J1,2)+UtyS*DS(I1,J1,5)+UnS*DS(I1,J1,8)         
!         VZS=UtxS*DS(I1,J1,3)+UtyS*DS(I1,J1,6)+UnS*DS(I1,J1,9)         
!                                                                       
!         VXD=UtxD*DS(I1,J1,1)+UtyD*DS(I1,J1,4)+UnD*DS(I1,J1,7)         
!         VYD=UtxD*DS(I1,J1,2)+UtyD*DS(I1,J1,5)+UnD*DS(I1,J1,8)         
!         VZD=UtxD*DS(I1,J1,3)+UtyD*DS(I1,J1,6)+UnD*DS(I1,J1,9)         
!                                                                       
          VX =Utx*DS(I1,J1,1) + Uty*DS(I1,J1,4) + Un*DS(I1,J1,7) 
          VY =Utx*DS(I1,J1,2) + Uty*DS(I1,J1,5) + Un*DS(I1,J1,8) 
          VZ =Utx*DS(I1,J1,3) + Uty*DS(I1,J1,6) + Un*DS(I1,J1,9) 
!                                                                       
          U = U+VX 
          V = V+VY 
          W = W+VZ 
!                                                                       
               !J1                                                      
        ENDDO 
!                                                                       
               !I1                                                      
      ENDDO 
!                                                                       
      RETURN 
          !VELOCE2                                                      
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VEL_WINGFF(CRR,XC,YC,ZC,DUB,SIG,VX,VY,VZ,              &
     &                    VXS,VYS,VZS,VXD,VYD,VZD,PNS2)                 
! come VEL_WING ma prova FAR FIELD ctb                                  
!     CALCOLA velocita' sulla scia dovuta a sorgenti/doppiette sul corpo
!     da PROGRAM No. 12: INFLUENCE COEFF. OF A RECTILINEAR SOURCE/DOUBLE
!     BASED on sec 10.4.1 (pag. 245) and 10.4.2                         
!     ----------------------------------------------------------------- 
!     THIS PROGRAM CALCULATES THE INFLUENCE OF A RECTILINEAR PANEL AT AN
!     ARBITRARY POINT XC,YC,ZC                                          
      DIMENSION CRR(12) 
      DIMENSION X(5),Y(5),Z(5) 
      COMMON/NO3/ FOURPAY 
!                                                                       
!     EPS=1.E-12                                                        
!     EPS=1.E-06                                                        
      ONE=1. 
      EPS=EPSILON(ONE)*10. 
!                                                                       
!     PNDS=1.0                                                          
      SSIG=SIG/FOURPAY 
      DDUB=DUB/FOURPAY 
!                                                                       
!     SQUARE/FLAT PANEL                                                 
!                         !  4------------3                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  1------------2                             
      X(1) = CRR(1) 
      Y(1) = CRR(2) 
      Z(1) = CRR(3) 
      X(2) = CRR(4) 
      Y(2) = CRR(5) 
      Z(2) = CRR(6) 
      X(3) = CRR(7) 
      Y(3) = CRR(8) 
      Z(3) = CRR(9) 
      X(4) = CRR(10) 
      Y(4) = CRR(11) 
      Z(4) = CRR(12) 
      X(5)=X(1) 
      Y(5)=Y(1) 
      Z(5)=Z(1) 
!     CENTROIDE                                                         
!     PXO=(X(1)+X(2)+X(3)+X(4))/4.                                      
!     PYO=(Y(1)+Y(2)+Y(3)+Y(4))/4.                                      
!     PZO=(Z(1)+Z(2)+Z(3)+Z(4))/4.                                      
!                                                                       
      VXS = 0.0 
      VYS = 0.0 
      VZS = 0.0 
      VXD = 0.0 
      VYD = 0.0 
      VZD = 0.0 
!                                                                       
            !PXI-PXO                                                    
      PX=XC 
            !PYI-PYO                                                    
      PY=YC 
            !PZI-PZO                                                    
      PZ=ZC 
!                                                                       
!     RDIST=SQRT(PX*PX+PY*PY+PZ*PZ)                                     
!                                                                       
      PNLX=.25*(X(1)+X(2)+X(3)+X(4)) 
      PNLY=.25*(Y(1)+Y(2)+Y(3)+Y(4)) 
      PNLZ=.25*(Z(1)+Z(2)+Z(3)+Z(4)) 
! se tolgo calcolo far-field non servono PNS,PNX,PNY                    
      PNX=PX-PNLX 
      PNY=PY-PNLY 
      PNZ=PZ-PNLZ 
      PNS2=PNX*PNX+PNY*PNY+PNZ*PNZ 
!                                                                       
      D1X=X(3)-X(1) 
      D1Y=Y(3)-Y(1) 
      D1Z=Z(3)-Z(1) 
      D2X=X(4)-X(2) 
      D2Y=Y(4)-Y(2) 
      D2Z=Z(4)-Z(2) 
      CRX=D1Y*D2Z-D2Y*D1Z 
      CRY=D2X*D1Z-D1X*D2Z 
      CRZ=D1X*D2Y-D2X*D1Y 
      CRSQ=SQRT(CRX*CRX+CRY*CRY+CRZ*CRZ) 
                   !serve x calcolo ctb far-field                       
      AREA=CRSQ/2. 
! FAR FIELD OK x PNS2 > 9 (diff 10-7)                                   
      IF(PNS2.GT.9.) GO TO 100 
!                                                                       
      CNX=CRX/CRSQ 
      CNY=CRY/CRSQ 
      CNZ=CRZ/CRSQ 
!     PNN=CNX*PNX+CNY*PNY+CNZ*PNZ                                       
      TCMX=(X(3)+X(4))*.5 - PNLX 
      TCMY=(Y(3)+Y(4))*.5 - PNLY 
      TCMZ=(Z(3)+Z(4))*.5 - PNLZ 
      TMS=SQRT(TCMX*TCMX+TCMY*TCMY+TCMZ*TCMZ) 
!     CMX=((X(3)+X(4))/2. - PNLX)/TMS                                   
!     CMY=((Y(3)+Y(4))/2. - PNLY)/TMS                                   
!     CMZ=((Z(3)+Z(4))/2. - PNLZ)/TMS                                   
      CMX=TCMX/TMS 
      CMY=TCMY/TMS 
      CMZ=TCMZ/TMS 
      CLX=CMY*CNZ-CNY*CMZ 
      CLY=CNX*CMZ-CMX*CNZ 
      CLZ=CMX*CNY-CNX*CMY 
!                                                                       
!     WDNF=0.                                                           
      DO J=1,4 
         K=J+1 
         AX=PX-X(J) 
         AY=PY-Y(J) 
         AZ=PZ-Z(J) 
         BX=PX-X(K) 
         BY=PY-Y(K) 
         BZ=PZ-Z(K) 
         SX=X(K)-X(J) 
         SY=Y(K)-Y(J) 
         SZ=Z(K)-Z(J) 
         A=SQRT(AX*AX + AY*AY + AZ*AZ) 
         B=SQRT(BX*BX + BY*BY + BZ*BZ) 
         S=SQRT(SX*SX + SY*SY + SZ*SZ) 
!                                                                       
         IF(ABS(S).LT.EPS)write(60,*)                                   &
     &    'S <1E-5 suv VEL_WING'                                        
!                                                                       
!        SOURCE CONTRIBUTION                                            
!                                                                       
         SM=SX*CMX+SY*CMY+SZ*CMZ 
         SL=SX*CLX+SY*CLY+SZ*CLZ 
         AM=AX*CMX+AY*CMY+AZ*CMZ 
         AL=AX*CLX+AY*CLY+AZ*CLZ 
         ALL=AM*SL-AL*SM 
         IF((A+B-S).GT.0.0.AND.S.GT.0.0)THEN 
           RJ3=ALOG((A+B+S)/(A+B-S))/S 
         ELSE 
           RJ3=0.0 
         ENDIF 
         PA=PNZ*PNZ*SL + ALL*AM 
         PB=PA - ALL*SM 
         RNUM=SM*PNZ*(B*PA - A*PB) 
         DNOM=PA*PB + PNZ*PNZ*A*B*SM*SM 
         IF(ABS(PNZ).LT.EPS)THEN 
           DE=0.0 
         ELSE 
           IF(ABS(RNUM).GT.EPS)THEN 
             DE=ATAN2(RNUM,DNOM) 
           ELSE 
             DE=0.0 
           ENDIF 
         ENDIF 
!                                                                       
         VXS = VXS+(RJ3*(SM*CLX-SL*CMX)+DE*CNX) 
         VYS = VYS+(RJ3*(SM*CLY-SL*CMY)+DE*CNY) 
         VZS = VZS+(RJ3*(SM*CLZ-SL*CMZ)+DE*CNZ) 
!                                                                       
!        DOUBLET CONTRIBUTION                                           
!                                                                       
                              !A vec B!_x                               
         AVBX = AY*BZ - AZ*BY 
                              !A vec B!_y                               
         AVBY = AZ*BX - AX*BZ 
                              !A vec B!_z                               
         AVBZ = AX*BY - AY*BX 
                                     !A dot B                           
         ADB = AX*BX + AY*BY + AZ*BZ 
         AB=A*B 
         if(abs(AB).LE.0.000001) stop 'A*B =0 sub VEL_WING' 
         VMOD=(A+B)/(AB*(AB + ADB)) 
         VXD = VXD + VMOD*AVBX 
         VYD = VYD + VMOD*AVBY 
         VZD = VZD + VMOD*AVBZ 
!        NEAR FIELD DOUBLET w-CONTRIBUTION WDNF (U=V=0)                 
!        RJ=SQRT(AX*AX+AY*AY)                                           
!        RK=SQRT(BX*BX+BY*BY) !J+1                                      
!        SX=X(K)-X(J)                                                   
!        SY=Y(K)-Y(J)                                                   
!        WDNF= WDNF + ((RJ+RK)*(AY*BX-AX*BY))/                          
!    &   (RJ*RK*(RJ*RK-(AX*BX+AY*BY)))                                  
            !J                                                          
      ENDDO 
!                                                                       
      VXS =  SSIG*VXS 
      VYS =  SSIG*VYS 
      VZS =  SSIG*VZS 
      VXD =  DDUB*VXD 
      VYD =  DDUB*VYD 
      VZD =  DDUB*VZD 
!     WDNF = DDUB*WDNF                                                  
!                                                                       
      VX = VXD + VXS 
      VY = VYD + VYS 
      VZ = VZD + VZS 
      RETURN 
!                                                                       
! calcolo FAR-FIELD ctb                                                 
! x FAR-FIELD CONTRIBUTIONS see (10.111-10.113) x DOUBLET ecc.          
  100 PNS=SQRT(PNS2) 
      PNS3=PNS2*PNS 
      PNS5=PNS3*PNS2 
      PHIS_FAR=-SSIG*AREA/PNS 
      TCM=-PHIS_FAR/PNS2 
      VXS=TCM*PNX 
      VYS=TCM*PNY 
      VZS=TCM*PNZ 
!                                                                       
      TCM=DDUB*AREA/PNS5 
      TCM2=3.*TCM*PNZ 
      VXD=TCM2*PNX 
      VYD=TCM2*PNY 
!     VZD=-TCM*(PNX*PNX+PNY*PNY-2.*PNZ*PNZ)                             
      VZD=-TCM*(PNS2-3.*PNZ*PNZ) 
!                                                                       
      VX = VXD + VXS 
      VY = VYD + VYS 
      VZ = VZD + VZS 
!                                                                       
      RETURN 
             !VEL_WINGFF                                                
      END                                           
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VEL_WAKE(JB,NSTEPS,X,Y,Z,QFW,DUBW,N,U,V,W) 
!                                                                       
!     CALCULATES VELOCITY AT A COLLOCATION POINT X,Y,Z                  
!     DUE TO WAKE DOUBLET DISTRIBUTION DUBW(N,JJ) IN A WING FIXED       
!     COORDINATE SYSTEM AT T=N*DT                                       
!                                                                       
                          ! coordin. centro generico pannello scia      
      DIMENSION QCW(3) 
                                      !corner points inizio pannelli sci
      DIMENSION QFW(NSTEPS+1,JB+1,3) 
                                !intensità doppiette scia              
      DIMENSION DUBW(NSTEPS,JB) 
                          ! vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                          !normale e tangenti a pannello generico corpo 
      DIMENSION DSS(10) 
!                                                                       
!     X,Y,Z = WAKE COLLOCATION POINTS                                   
      U=0. 
      V=0. 
      W=0. 
      DO IT=1,N 
        DO JJ=1,JB 
!       IT,JJ = IFLUENCING PANEL in WAKE                                
          DO K=1,3 
            QCW(K)=(QFW(IT+1,JJ,K)+QFW(IT+1,JJ+1,K)+                    &
     &              QFW(IT,JJ,K)+QFW(IT,JJ+1,K))/4.                     
                !K                                                      
          ENDDO 
!                                                                       
!     COMPUTATION OF TANGENTIAL CHORDWISE AND SPANWISE VECTORS          
!     DSS(1-3) and DSS(4-6)                                             
!     NORMAL VECTOR DSS(7-9), PANEL AREA DSS(10)                        
!     calcola il DSS del pannello IT,JJ                                 
!     le coordin della scia QFW sono nel SR corpo e anche il DSS        
!                                                                       
          CALL PANEL(                                                   &
     &    QFW(IT+1,JJ,1),QFW(IT+1,JJ,2),QFW(IT+1,JJ,3),                 &
     &    QFW(IT,JJ,1),QFW(IT,JJ,2),QFW(IT,JJ,3),                       &
     &    QFW(IT+1,JJ+1,1),QFW(IT+1,JJ+1,2),QFW(IT+1,JJ+1,3),           &
     &    QFW(IT,JJ+1,1),QFW(IT,JJ+1,2),QFW(IT,JJ+1,3),                 &
     &    DSS(1),DSS(2),DSS(3),DSS(4),DSS(5),                           &
     &    DSS(6),DSS(7),DSS(8),DSS(9),DSS(10))                          
!                                                                       
!         CONVERT COLLOCATION POINT TO PANEL COORDINATES                
          XB=X-QCW(1) 
          YB=Y-QCW(2) 
          ZB=Z-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,XC,YC,ZC) 
!                                                                       
!         TRANSFORMATION OF A WAKE POINT QCW(1),QCW(2),QCW(3)           
!         INTO PANEL COORDINATES CRR(1),CRR(2),CRR(3)                   
!         QFW(I,J,1),QFW(I,J,2),QFW(I,J,3) = WAKE PANEL CORNERPOINT     
!         CRR = coordinate nel SR ala con centro in QCW                 
!                                                                       
          XB=QFW(IT+1,JJ,1)-QCW(1) 
          YB=QFW(IT+1,JJ,2)-QCW(2) 
          ZB=QFW(IT+1,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(1),CRR(2),CRR(3)) 
!                                                                       
          XB=QFW(IT,JJ,1)-QCW(1) 
          YB=QFW(IT,JJ,2)-QCW(2) 
          ZB=QFW(IT,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(4),CRR(5),CRR(6)) 
!                                                                       
          XB=QFW(IT,JJ+1,1)-QCW(1) 
          YB=QFW(IT,JJ+1,2)-QCW(2) 
          ZB=QFW(IT,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(7),CRR(8),CRR(9)) 
!                                                                       
          XB=QFW(IT+1,JJ+1,1)-QCW(1) 
          YB=QFW(IT+1,JJ+1,2)-QCW(2) 
          ZB=QFW(IT+1,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(10),CRR(11),CRR(12)) 
!                                                                       
!         CALL P12VEL(CRR,XC,YC,ZC,DUBW(IT,JJ),Utx,Uty,Un)              
! alternativa con calcolo Far-Field                                     
          CALL VEL_DOUBFF(CRR,XC,YC,ZC,DUBW(IT,JJ),Utx,Uty,Un,PNS2) 
!c   &                    UtxFF,UtyFF,UnFF,PNS2)                        
!                                                                       
!         write(137,*)PNS2,Utx-UtxFF,Uty-UtyFF,Un-UnFF,                 
!    &               'PNS2,Utx-UtxFF,Uty-UtyFF,Un-UnFF,'                
!                                                                       
!         converte da SR pannelo a SR corpo                             
          VX =Utx*DSS(1) + Uty*DSS(4) + Un*DSS(7) 
          VY =Utx*DSS(2) + Uty*DSS(5) + Un*DSS(8) 
          VZ =Utx*DSS(3) + Uty*DSS(6) + Un*DSS(9) 
!                                                                       
          U = U + VX 
          V = V + VY 
          W = W + VZ 
!                                                                       
               !JJ                                                      
        ENDDO 
               !IT                                                      
      ENDDO 
!                                                                       
      RETURN 
          !VEL_WAKE                                                     
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE INFLUENCEW(DOUB,XC,YC,ZC,CRR) 
      DIMENSION CRR(12) 
      COMMON/NO3/ FOURPAY 
!                                                                       
!     DOUBLET POTENTIAL AT POINT (XC,YC,ZC) DUE TO                      
!     PANEL I,J  i.e. CR(I,J,1-2-3) (X1,Y1,Z1,...X4,Y4,Z4),             
!     SEE KATZ & PLOTKIN PP 248 (10.105)                                
!     ?? SEE KATZ & PLOTKIN PP 283-6                                    
!     SEE ALSO Hess & Smith 1967 p 50-58                                
!     SEE ALSO Ortega et al. 2010 p 11                                  
!     SEE ALSO Maskew 1987 p 37                                         
!     risultati = a sub P12W                                            
!                                                                       
!     EP=0.0000000001                                                   
      ONE=1. 
      EP=EPSILON(ONE)*10. 
!                                                                       
      DOUB=0.0 
!     IF(ABS(ZC).LT.EP) THEN    !rivedere se serve                      
!C      DOUB=-0.5                                                       
!       write(28,*)'ABS(ZC).LT.EP in sub INFLUENCEW'                    
!     ENDIF                                                             
!                                                                       
      X1 = CRR(1) 
      Y1 = CRR(2) 
      Z1 = CRR(3) 
      X2 = CRR(4) 
      Y2 = CRR(5) 
      Z2 = CRR(6) 
      X3 = CRR(7) 
      Y3 = CRR(8) 
      Z3 = CRR(9) 
      X4 = CRR(10) 
      Y4 = CRR(11) 
      Z4 = CRR(12) 
!                                                                       
!     EP, PANEL SIDE CUTOFF DISTANCE                                    
!     PANEL SIDE (D) DISTANCE (R), E, AND H (EQS. 10.90 & 10.92-10.94)  
!                                                                       
      X2mX1=X2-X1 
      X3mX2=X3-X2 
      X4mX3=X4-X3 
      X1mX4=X1-X4 
      Y2mY1=Y2-Y1 
      Y3mY2=Y3-Y2 
      Y4mY3=Y4-Y3 
      Y1mY4=Y1-Y4 
      XCmX1=XC-X1 
      XCmX2=XC-X2 
      XCmX3=XC-X3 
      XCmX4=XC-X4 
      YCmY1=YC-Y1 
      YCmY2=YC-Y2 
      YCmY3=YC-Y3 
      YCmY4=YC-Y4 
!                                                                       
                                      !d12 (10.90a)                     
      D1=SQRT((X2mX1)**2+(Y2mY1)**2) 
                                      !d23 (10.90a)                     
      D2=SQRT((X3mX2)**2+(Y3mY2)**2) 
                                      !d34 (10.90c)                     
      D3=SQRT((X4mX3)**2+(Y4mY3)**2) 
                                      !d41 (10.90d)                     
      D4=SQRT((X1mX4)**2+(Y1mY4)**2) 
!                                                                       
                                             !(10.92)                   
      R1=SQRT((XCmX1)**2+(YCmY1)**2+ZC**2) 
                                             !(10.92)                   
      R2=SQRT((XCmX2)**2+(YCmY2)**2+ZC**2) 
                                             !(10.92)                   
      R3=SQRT((XCmX3)**2+(YCmY3)**2+ZC**2) 
                                             !(10.92)                   
      R4=SQRT((XCmX4)**2+(YCmY4)**2+ZC**2) 
!                                                                       
                           !(10.93)                                     
      E1=(XCmX1)**2+ZC**2 
                           !(10.93)                                     
      E2=(XCmX2)**2+ZC**2 
                           !(10.93)                                     
      E3=(XCmX3)**2+ZC**2 
                           !(10.93)                                     
      E4=(XCmX4)**2+ZC**2 
!                                                                       
                           !(10.94)                                     
      H1=(XCmX1)*(YCmY1) 
                           !(10.94)                                     
      H2=(XCmX2)*(YCmY2) 
                           !(10.94)                                     
      H3=(XCmX3)*(YCmY3) 
                           !(10.94)                                     
      H4=(XCmX4)*(YCmY4) 
!                                                                       
!     DOUBLET (Q, DOUB) INFLUENCE IN PANEL COORDINATES                  
!     FOR TRIANGULAR PANEL THE 4TH SIDE CONTRIBUTION IS ZERO            
!                                                                       
      Q1=0. 
      IF (D1.GE.EP) THEN 
                                !PA                                     
        F=(Y2mY1)*E1-(X2mX1)*H1 
                                !PB                                     
        G=(Y2mY1)*E2-(X2mX1)*H2 
!       R2=B; R1=A; ZC=PN; X2-X1=SM; D1=S                               
        RNUM=ZC*(X2mX1)*(F*R2-G*R1) 
        DNOM=(ZC*X2mX1)**2*R1*R2 + F*G 
                            !(10.105)                                   
        Q1=ATAN2(RNUM,DNOM) 
      ENDIF 
!                                                                       
      Q2=0. 
      IF (D2.GE.EP) THEN 
        F=(Y3mY2)*E2-(X3mX2)*H2 
        G=(Y3mY2)*E3-(X3mX2)*H3 
        RNUM=ZC*(X3mX2)*(F*R3-G*R2) 
        DNOM=(ZC*X3mX2)**2*R2*R3 + F*G 
                            !(10.105)                                   
        Q2=ATAN2(RNUM,DNOM) 
      ENDIF 
!                                                                       
      Q3=0. 
      IF (D3.GE.EP) THEN 
        F=(Y4mY3)*E3-(X4mX3)*H3 
        G=(Y4mY3)*E4-(X4mX3)*H4 
        RNUM=ZC*(X4mX3)*(F*R4-G*R3) 
        DNOM=(ZC*X4mX3)**2*R3*R4 + F*G 
                            !(10.105)                                   
        Q3=ATAN2(RNUM,DNOM) 
      ENDIF 
!                                                                       
      Q4=0. 
      IF (D4.GE.EP) THEN 
        F=(Y1mY4)*E4-(X1mX4)*H4 
        G=(Y1mY4)*E1-(X1mX4)*H1 
        RNUM=ZC*(X1mX4)*(F*R1-G*R4) 
        DNOM=(ZC*X1mX4)**2*R4*R1 + F*G 
                            !(10.105)                                   
        Q4=ATAN2(RNUM,DNOM) 
      ENDIF 
!                                                                       
!     ADD CONTRIBUTIONS FROM THE 4 SIDES                                
!                                                                       
                                             !times doublet strength (10
      DOUB=+(Q1+Q2+Q3+Q4)/FOURPAY 
!                                                                       
      RETURN 
          !INFLUENCEW                                                   
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VELOCE3(IB,JB,QC,CR,DUB,SIGMA,DS,XX,YY,ZZ,U,V,W) 
!                                                                       
!     CALCULATES WAKE VELOCITY AT A COLLOCATION POINT XX,YY,ZZ i.e. QFW(
!     DUE TO WING SOURCE AND DOUBLET DISTRIBUTION ON PANEL (I1,J1)      
!     IN A WING FIXED COORDINATE SYSTEM                                 
!                                                                       
                               !coordinate centro pannelli corpo        
      DIMENSION QC(IB+1,JB,3) 
                               !coordinate 4 vertici pannelli corpo     
      DIMENSION CR(IB+1,JB,12) 
                             !doppiette su corpo                        
      DIMENSION DUB(IB,JB) 
                             !sorgenti su corpo                         
      DIMENSION SIGMA(IB,JB) 
                               !normale e tangenti a pannello corpo     
      DIMENSION DS(IB+1,JB,10) 
                           !vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                           !normale e tangenti a pannello generico corpo
      DIMENSION DSS(10) 
!                                                                       
      COMMON/NO3/ FOURPAY 
!                                                                       
!     WING FIXED PANEL CORNERPOINTS                                     
!     CR(I,J,1)=X(i,j); CR(I,J,2)=Y(I,J); CR(I,J,3)=Z(I,J)              
! calcola il potenziale in P(x,y,z) prodotta dal pannello I,J           
!     i,j I---1----I i,j+1                                              
!         I        I                                                    
!         4  I,J   2                         P(x,y,z)                   
!         I        I                                                    
!   i+1,j I---3----I i+1,j+1                                            
!                                                                       
      IB1=IB+1 
!                                                                       
!-riv X1=(X-X0)*CS1+(Z-Z0)*SN1                                          
!-riv Y1=Y                                                              
!-riv Z1=-(X-X0)*SN1+(Z-Z0)*CS1                                         
!     X,Y,Z è il punto della scia ove voglio calcolare la velocità    
!                                                                       
      U=0. 
      V=0. 
      W=0. 
!                                                                       
!     INFLUENCING PANEL I1,J1 sul corpo                                 
      DO I1=1,IB 
        DO J1=1,JB 
!         err1=0.                                                       
!                                                                       
          DO K=1,10 
            DSS(K)=DS(I1,J1,K) 
          ENDDO 
!         CONVERT COLLOCATION POINT XX,YY,ZZ TO PANEL (I1,J1) COORDINATE
!  primi 3= coordin centro influencing panel (I1,J1),                   
!  secondi 3= colloc pt (I,J) da trasformare; X,Y,Z= nodi della scia QFW
!  XC,YC,ZC= colloc pt nel sr pannello                                  
          X = XX - QC(I1,J1,1) 
          Y = YY - QC(I1,J1,2) 
          Z = ZZ - QC(I1,J1,3) 
                                                                        
          CALL CONVERT1(X,Y,Z,DSS,XC,YC,ZC) 
                                                                        
!         effetto del pannello I1,J1 sul collocation pt X,Y,Z           
          DO K=1,12 
            CRR(K)=CR(I1,J1,K) 
          ENDDO 
!                                                                       
!         CALL VEL_WING(CRR,XC,YC,ZC,DUB(I1,J1),SIGMA(I1,J1),           
!    &                  Utx,Uty,Un,UtxS,UtyS,UnS,UtxD,UtyD,UnD)         
          CALL VEL_WINGFF(CRR,XC,YC,ZC,DUB(I1,J1),SIGMA(I1,J1),         &
     &    FUtx,FUty,FUn,FUtxS,FUtyS,FUnS,FUtxD,FUtyD,FUnD,PNS2)         
!                                                                       
!         err1=max(ERR1,ABS(Utx-FUtx))                                  
!         err1=max(ERR1,ABS(Uty-FUty))                                  
!         err1=max(ERR1,ABS(Un-FUn))                                    
!        write(136,*) I1,J1,PNS2,ERR1,'INFLUENCING PANEL I1,J1 su corpo'
!         write(136,*)Utx-FUtx,Uty-FUty,Un-FUn,                         
!    &               'Utx-FUtx,Uty-FUty,Un-FUn'                         
!         write(136,*)UtxS-FUtxS,UtyS-FUtyS,UnS-FUnS,                   
!    &               'UtxS-FUtxS,UtyS-FUtyS,UnS-FUnS'                   
!         write(136,*)UtxD-FUtxD,UtyD-FUtyD,UnD-FUnD,                   
!    &               'UtxD-FUtxD,UtyD-FUtyD,UnD-FUnD'                   
!                                                                       
!         proietta per calcolare le componenti nel S.R corpo            
!         ux=utx*tx_x+uty*ty_x+un*n_x                                   
!         uy=utx*tx_y+uty*ty_y+un*n_y                                   
!         uz=utx*tx_z+uty*ty_z+un*n_z                                   
!                                                                       
!         VXS=UtxS*DS(I1,J1,1)+UtyS*DS(I1,J1,4)+UnS*DS(I1,J1,7)         
!         VYS=UtxS*DS(I1,J1,2)+UtyS*DS(I1,J1,5)+UnS*DS(I1,J1,8)         
!         VZS=UtxS*DS(I1,J1,3)+UtyS*DS(I1,J1,6)+UnS*DS(I1,J1,9)         
!                                                                       
!         VXD=UtxD*DS(I1,J1,1)+UtyD*DS(I1,J1,4)+UnD*DS(I1,J1,7)         
!         VYD=UtxD*DS(I1,J1,2)+UtyD*DS(I1,J1,5)+UnD*DS(I1,J1,8)         
!         VZD=UtxD*DS(I1,J1,3)+UtyD*DS(I1,J1,6)+UnD*DS(I1,J1,9)         
!                                                                       
          VX =Utx*DS(I1,J1,1) + Uty*DS(I1,J1,4) + Un*DS(I1,J1,7) 
          VY =Utx*DS(I1,J1,2) + Uty*DS(I1,J1,5) + Un*DS(I1,J1,8) 
          VZ =Utx*DS(I1,J1,3) + Uty*DS(I1,J1,6) + Un*DS(I1,J1,9) 
!                                                                       
          U = U+VX 
          V = V+VY 
          W = W+VZ 
!                                                                       
! per cfr con 2D devo guardare i ctb di tutta la riga J                 
!                                                                       
               !J1                                                      
        ENDDO 
!                                                                       
               !I1                                                      
      ENDDO 
!     write(136,*) ERR1,'ERR1'                                          
!                                                                       
      RETURN 
          !VELOCE3                                                      
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VEL_WAKE3(JB,NSTEPS,X,Y,Z,QFW,DUBW,N,U,V,W) 
! come VEL_WAKE ma con scritte x check                                  
!     CALCULATES VELOCITY AT A COLLOCATION POINT X,Y,Z                  
!     DUE TO WAKE DOUBLET DISTRIBUTION DUBW(N,JJ) IN A WING FIXED       
!     COORDINATE SYSTEM AT T=N*DT                                       
!                                                                       
                          ! coordin. centro generico pannello scia      
      DIMENSION QCW(3) 
                                      !corner points inizio pannelli sci
      DIMENSION QFW(NSTEPS+1,JB+1,3) 
                                !intensità doppiette scia              
      DIMENSION DUBW(NSTEPS,JB) 
                          ! vettore contenente corner points di un panne
      DIMENSION CRR(12) 
                          !normale e tangenti a pannello generico corpo 
      DIMENSION DSS(10) 
!                                                                       
!     X,Y,Z = WAKE COLLOCATION POINTS                                   
      U=0. 
      V=0. 
      W=0. 
          DIFMX=0. 
      DO IT=1,N 
        DO JJ=1,JB 
!       IT,JJ = IFLUENCING PANEL in WAKE                                
          DO K=1,3 
            QCW(K)=(QFW(IT+1,JJ,K)+QFW(IT+1,JJ+1,K)+                    &
     &              QFW(IT,JJ,K)+QFW(IT,JJ+1,K))/4.                     
                !K                                                      
          ENDDO 
!                                                                       
!     COMPUTATION OF TANGENTIAL CHORDWISE AND SPANWISE VECTORS          
!     DSS(1-3) and DSS(4-6)                                             
!     NORMAL VECTOR DSS(7-9), PANEL AREA DSS(10)                        
!     calcola il DSS del pannello IT,JJ                                 
!                                                                       
          CALL PANEL(                                                   &
     &    QFW(IT+1,JJ,1),QFW(IT+1,JJ,2),QFW(IT+1,JJ,3),                 &
     &    QFW(IT,JJ,1),QFW(IT,JJ,2),QFW(IT,JJ,3),                       &
     &    QFW(IT+1,JJ+1,1),QFW(IT+1,JJ+1,2),QFW(IT+1,JJ+1,3),           &
     &    QFW(IT,JJ+1,1),QFW(IT,JJ+1,2),QFW(IT,JJ+1,3),                 &
     &    DSS(1),DSS(2),DSS(3),DSS(4),DSS(5),                           &
     &    DSS(6),DSS(7),DSS(8),DSS(9),DSS(10))                          
!                                                                       
!         CONVERT COLLOCATION POINT TO PANEL COORDINATES                
          XB=X-QCW(1) 
          YB=Y-QCW(2) 
          ZB=Z-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,XC,YC,ZC) 
!                                                                       
!         TRANSFORMATION OF A WAKE POINT QCW(1),QCW(2),QCW(3)           
!         INTO PANEL COORDINATES CRR(1),CRR(2),CRR(3)                   
!         QFW(I,J,1),QFW(I,J,2),QFW(I,J,3) = WAKE PANEL CORNERPOINT     
!         CRR = coordinate nel SR ala con centro in QCW                 
!                                                                       
          XB=QFW(IT+1,JJ,1)-QCW(1) 
          YB=QFW(IT+1,JJ,2)-QCW(2) 
          ZB=QFW(IT+1,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(1),CRR(2),CRR(3)) 
!                                                                       
          XB=QFW(IT,JJ,1)-QCW(1) 
          YB=QFW(IT,JJ,2)-QCW(2) 
          ZB=QFW(IT,JJ,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(4),CRR(5),CRR(6)) 
!                                                                       
          XB=QFW(IT,JJ+1,1)-QCW(1) 
          YB=QFW(IT,JJ+1,2)-QCW(2) 
          ZB=QFW(IT,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(7),CRR(8),CRR(9)) 
!                                                                       
          XB=QFW(IT+1,JJ+1,1)-QCW(1) 
          YB=QFW(IT+1,JJ+1,2)-QCW(2) 
          ZB=QFW(IT+1,JJ+1,3)-QCW(3) 
          CALL CONVERT1(XB,YB,ZB,DSS,CRR(10),CRR(11),CRR(12)) 
!                                                                       
!         subr. alternativa sostituita da VEL_DOUBFF cha ha calcolo farf
!         CALL P12VEL(CRR,XC,YC,ZC,DUBW(IT,JJ),Utx,Uty,Un)              
!                                                                       
          CALL VEL_DOUBFF(CRR,XC,YC,ZC,DUBW(IT,JJ),                     &
     &                    UtxFF,UtyFF,UnFF,PNS2)                        
!                                                                       
          DIFMX=MAX(DIFMX,ABS(Utx-UtxFF)) 
          DIFMX=MAX(DIFMX,ABS(Uty-UtyFF)) 
          DIFMX=MAX(DIFMX,ABS(Un-UnFF)) 
!         write(137,*)PNS2,Utx-UtxFF,Uty-UtyFF,Un-UnFF,                 
!    &               'PNS2,Utx-UtxFF,Uty-UtyFF,Un-UnFF,'                
!         write(137,*)IT,JJ,PNS2,Uty-UtyFF,Uty,UtyFF,                   
!    &               'IT,JJ,PNS2,Uty-UtyFF,Uty,UtyFF'                   
!         write(137,*)IT,JJ,PNS2,Un-UnFF,Un,UnFF,                       
!    &               'IT,JJ,PNS2,Un-UnFF,Un,UnFF'                       
!                                                                       
!         converte da SR pannelo a SR corpo                             
          VX =Utx*DSS(1) + Uty*DSS(4) + Un*DSS(7) 
          VY =Utx*DSS(2) + Uty*DSS(5) + Un*DSS(8) 
          VZ =Utx*DSS(3) + Uty*DSS(6) + Un*DSS(9) 
!                                                                       
          U = U + VX 
          V = V + VY 
          W = W + VZ 
!                                                                       
               !JJ                                                      
        ENDDO 
               !IT                                                      
      ENDDO 
      write(137,*)DIFMX,'DIFMX' 
!                                                                       
      RETURN 
          !VEL_WAKE3                                                    
      END                                           
!                                                                       
!                                                                       
!--------+---------+---------+---------+---------+---------+---------+--
      SUBROUTINE VEL_DOUBFF(CRR,XC,YC,ZC,DUB,VX,VY,VZ,PNS2) 
! come VEL_WINGFF ma con solo ctb doppiette x calcolo Vel scia su scia  
!     CALCOLA velocita' sulla scia dovuta a doppiette della scia        
!     BASED on sec 10.4.1 (pag. 245) and 10.4.2                         
!     ----------------------------------------------------------------- 
!     THIS PROGRAM CALCULATES THE INFLUENCE OF A RECTILINEAR PANEL AT AN
!     ARBITRARY POINT XC,YC,ZC                                          
      DIMENSION CRR(12) 
      DIMENSION X(5),Y(5),Z(5) 
      COMMON/NO3/ FOURPAY 
!                                                                       
!     EPS=1.E-12                                                        
!     EPS=1.E-06                                                        
      ONE=1. 
      EPS=EPSILON(ONE)*10. 
!                                                                       
!     PNDS=1.0                                                          
      DDUB=DUB/FOURPAY 
!                                                                       
!     SQUARE/FLAT PANEL                                                 
!                         !  4------------3                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  I            I                             
!                         !  1------------2                             
      X(1) = CRR(1) 
      Y(1) = CRR(2) 
      Z(1) = CRR(3) 
      X(2) = CRR(4) 
      Y(2) = CRR(5) 
      Z(2) = CRR(6) 
      X(3) = CRR(7) 
      Y(3) = CRR(8) 
      Z(3) = CRR(9) 
      X(4) = CRR(10) 
      Y(4) = CRR(11) 
      Z(4) = CRR(12) 
      X(5)=X(1) 
      Y(5)=Y(1) 
      Z(5)=Z(1) 
!     CENTROIDE                                                         
!     PXO=(X(1)+X(2)+X(3)+X(4))/4.                                      
!     PYO=(Y(1)+Y(2)+Y(3)+Y(4))/4.                                      
!     PZO=(Z(1)+Z(2)+Z(3)+Z(4))/4.                                      
!                                                                       
      VX = 0.0 
      VY = 0.0 
      VZ = 0.0 
!                                                                       
            !PXI-PXO                                                    
      PX=XC 
            !PYI-PYO                                                    
      PY=YC 
            !PZI-PZO                                                    
      PZ=ZC 
!                                                                       
!     RDIST=SQRT(PX*PX+PY*PY+PZ*PZ)                                     
!                                                                       
      PNLX=.25*(X(1)+X(2)+X(3)+X(4)) 
      PNLY=.25*(Y(1)+Y(2)+Y(3)+Y(4)) 
      PNLZ=.25*(Z(1)+Z(2)+Z(3)+Z(4)) 
      PNX=PX-PNLX 
      PNY=PY-PNLY 
      PNZ=PZ-PNLZ 
      PNS2=PNX*PNX+PNY*PNY+PNZ*PNZ 
! FAR FIELD OK x PNS2 > 9 (diff 10-7)                                   
      IF(PNS2.GT.9.) GO TO 100 
!                                                                       
      DO J=1,4 
         K=J+1 
         AX=PX-X(J) 
         AY=PY-Y(J) 
         AZ=PZ-Z(J) 
         BX=PX-X(K) 
         BY=PY-Y(K) 
         BZ=PZ-Z(K) 
                                       !R1                              
         A=SQRT(AX*AX + AY*AY + AZ*AZ) 
                                       !R2                              
         B=SQRT(BX*BX + BY*BY + BZ*BZ) 
!                                                                       
!        DOUBLET CONTRIBUTION                                           
!                                                                       
                              !(r1 \vec r2)_x                           
         AVBX = AY*BZ - AZ*BY 
                              !(r1 \vec r2)_y                           
         AVBY = AZ*BX - AX*BZ 
                              !(r1 \vec r2)_z                           
         AVBZ = AX*BY - AY*BX 
                                     !R1 dot R2                         
         ADB = AX*BX + AY*BY + AZ*BZ 
!                                                                       
         VMOD=0. 
         AB=A*B 
         DEN=AB*(AB + ADB) 
         if(abs(DEN).GT.EPS) VMOD=(A+B)/DEN 
!                                                                       
         VX = VX + VMOD*AVBX 
         VY = VY + VMOD*AVBY 
         VZ = VZ + VMOD*AVBZ 
            !J                                                          
      ENDDO 
      VX = VX *DDUB 
      VY = VY *DDUB 
      VZ = VZ *DDUB 
!                                                                       
      RETURN 
!                                                                       
! calcolo FAR-FIELD ctb                                                 
! x FAR-FIELD CONTRIBUTIONS see (10.111-10.113) x DOUBLET ecc.          
  100 CONTINUE 
      D1X=X(3)-X(1) 
!     D1Y=Y(3)-Y(1)                                                     
!     D1Z=Z(3)-Z(1)                                                     
!     D2X=X(4)-X(2)                                                     
      D2Y=Y(4)-Y(2) 
!     D2Z=Z(4)-Z(2)                                                     
!     CRX=D1Y*D2Z-D2Y*D1Z                                               
!     CRY=D2X*D1Z-D1X*D2Z                                               
!     CRZ=D1X*D2Y-D2X*D1Y                                               
!     CRSQ=SQRT(CRX*CRX+CRY*CRY+CRZ*CRZ)                                
!     AREA=CRSQ/2. !serve x calcolo ctb far-field                       
      AREA=D1X*D2Y 
!                                                                       
      PNS=SQRT(PNS2) 
      PNS5=PNS2*PNS2*PNS 
      TCM=DDUB*AREA/PNS5 
      TCM2=3.*TCM*PNZ 
      VX=TCM2*PNX 
      VY=TCM2*PNY 
!     VZ=-TCM*(PNX*PNX+PNY*PNY-2.*PNZ*PNZ)                              
      VZ=-TCM*(PNS2-3.*PNZ*PNZ) 
!                                                                       
      RETURN 
             !VEL_DOUBFF                                                
      END
!-----------------
      SUBROUTINE outfile(var,string)
character(len = 3)::string
integer:: var
!
if (var.le.9) then 
    write(string,'(i1)') abs(var)
else if (var.le.99) then
    write(string,'(i2)') var
else if (var.le.999) then
    write(string,'(i3)') var
end if
!
return
end subroutine
! ============================================================
SUBROUTINE CALC_PHIS(I,J,CR,PHIBAR)
real:: CR(:,:,:)
! formula valida solo per pannelli rettangolari planari
x1 = CR(I,J,1)
y1 = CR(I,J,2)
x2 = CR(I,J,4)
y2 = CR(I,J,5)
x3 = CR(I,J,7)
y3 = CR(I,J,8)
x4 = CR(I,J,10)
y4 = CR(I,J,11)
!
r = sqrt(x1**2 + y1**2)
!
phibar = x1*log((r + y2)/(r - y2)) +&
	 y2*log((r + x3)/(r - x3)) +&
	 x3*log((r + y4)/(r - y4)) +&
	 y4*log((r + x1)/(r - x1)) 
END SUBROUTINE CALC_PHIS
end module
