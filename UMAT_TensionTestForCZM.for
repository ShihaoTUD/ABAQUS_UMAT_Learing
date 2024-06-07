      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C     
      INCLUDE 'ABA_PARAM.INC'
C     
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),
     2     DDSDDT(NTENS),DRPLDE(NTENS),
     3     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0)
      DIMENSION STRANT(3)
      
      SIGMA_N=PROPS(1)
      SIGMA_S=PROPS(2)
      SIGMA_T=PROPS(3)

      E_N=PROPS(4)
      E_S=PROPS(5)
      E_T=PROPS(6)

      GIC=PROPS(7)
      GIIC=PROPS(8)
      GIIIC=PROPS(9)
      ETA=5.42D0
      IF(TIME(2).EQ.ZERO)THEN
          STATEV(1)=ZERO
      ENDIF

      DO I=1,3
          STRANT(I)=STRAN(I)+DSTRAN(I)
      ENDDO
      STRESS(1)=E_N*STRANT(1)
      STRESS(2)=E_S*STRANT(2)
      STRESS(3)=E_T*STRANT(3)
      
      DO I=1,NTENS
          DO J=1,NTENS
              DDSDDE(I,J)=ZERO
          ENDDO
      ENDDO
      DDSDDE(1,1)=E_N
      DDSDDE(2,2)=E_S
      DDSDDE(3,3)=E_T

C     BK LAW
      GC=GIC+(GIIC-GIC)*(GIIC/GIIIC)**ETA
C     ? the displacement and the strain ?
      DELTANO=SIGMA_N/E_N
      DELTANF=TWO*GIC/SIGMA_N
      DAMAGE=STATEV(1)
          IF(STRESS(1).GE.ZERO)THEN
              IF(STRANT(1).GE.ZERO.AND.STRANT(1).LE.DELTANO)THEN
                  DAMAGE=ZERO
                  STATEV(1)=MAX(STATEV(1),DAMAGE)
              ELSEIF(STRANT(1).GT.DELTANO.AND.STRANT(1).LT.DELTANF)THEN
                  DAMAGE=DELTANF*(STRANT(1)-DELTANO)/
     1                    (STRANT(1)*(DELTANF-DELTANO))
                  STATEV(1)=MAX(STATEV(1),DAMAGE)
              ELSE
                  STATEV(1)=ONE
              ENDIF
          ELSE
              STATEV(1)=ZERO
          ENDIF
              
      
      IF(STRANT(1).GT.DELTANO)THEN
          STRESS(1)=(ONE-STATEV(1))*E_N*STRANT(1)
          DDSDDE(1,1)=(ONE-STATEV(1))*E_N
      ELSE
          STRESS(1)=E_N*STRANT(1)
          DDSDDE(1,1)=E_N
      ENDIF
      
C   Controlling the element delection      
      IF(STATEV(1).GE.ONE)THEN
          STATEV(2)=ZERO
      ELSE
          STATEV(2)=ONE
      ENDIF
      STATEV(3)=STRESS(1)
      STATEV(4)=STRANT(1)
            
      RETURN
      END
      