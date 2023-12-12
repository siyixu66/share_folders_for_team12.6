      PROGRAM MoonPhaseCalcu
 !     	IMPLICIT NONE  ! use variable after definition

      !   Variable definition
              !For moon phase calculation
              implicit none
              INTEGER ::C  ! Central body
              INTEGER ::B  ! Target body
              REAL *8,DIMENSION(6) :: PO  !Coordinate vector:TMP
              REAL *8 :: t0   !the time when the Sun gives out light:TMP
              REAL *8 :: dt0,dt           !light time:TMP
              REAL *8,DIMENSION(3):: E2S  !Earth-Sun vector
              REAL *8,DIMENSION(3):: E2O  !Earth-Conical_point vector
              REAL *8,DIMENSION(3):: E2M  !Earth-Moon vector
              REAL *8,DIMENSION(3):: M2O  !Moon-Conical_point vector    
              REAL *8,DIMENSION(3):: S2M  !Sun-Moon vector   
              REAL *8 ::ANGLE_EO            !Earth cone Angle
              REAL *8 ::ANGLE_MO            !Moon cone Angle
              REAL *8 ::ANGLE_EOM           !Earth-Moon cone Angle
              REAL *8 ::theta               !Sun-Earth-Moon Angle:TMP
              REAL *8 ::PRODUCT             !Inner product
              REAL *8 ::ERR                 !Angular difference 
              REAL *8, PARAMETER::PI = 3.14159265358979323846
              REAL *8, PARAMETER::THRESHOLD = 175.0 * PI / 180.0    
              REAL *8::iteration_val=1.0/(24*3600)      !Iteration step(/day)
              REAL *8::light_time_diff=0.001/86400.0    !Threshold of Equation of light (/day)          
              
              !For time module
              !Gregorian Calendar_Start time
              INTEGER ::Year=2025            !CAL-year
              INTEGER ::Month=1              !CAL-month
              INTEGER ::Day=1                !CAL-day
              REAL *8 ::DDay=0               !CAL-dday(fractional part)    
              INTEGER ::Hour=0               !CAL-hour: UTC
              INTEGER ::Minute=0             !CAL-minute
              INTEGER ::Second=0             !CAL-second
              
              INTEGER ::NDP=3                        !number of decimal places of seconds
              CHARACTER :: SIGN                      !
              INTEGER ,DIMENSION(4) :: IHMSF         !hours, minutes, seconds, fraction
              !Gregorian Calendar_End time
              INTEGER ::Year_1=2026          !CAL-year
              INTEGER ::Month_1=1            !CAL-month
              INTEGER ::Day_1=1              !CAL-day
              REAL *8 ::DDay_1=0             !CAL-dday(fractional part)                  
              INTEGER ::Hour_1=0             !CAL-hour: UTC
              INTEGER ::Minute_1=0           !CAL-minute
              INTEGER ::Second_1=0           !CAL-second              
                
              !Julian Calendar:be in use in moon phase calculation
              !Julian Calendar_Start_time
              REAL *8 ::DJM0      !MJD zero-point: always 2400000.5D0
              REAL *8 ::DJM       !Modified Julian Date for 0 hrs
              !Julian Calendar_End_time
              REAL *8 ::DJM0_1    !MJD zero-point: always 2400000.5D0
              REAL *8 ::DJM_1     !Modified Julian Date for 0 hrs                        
              INTEGER ::J                  !State



              !Constants
              REAL *8 ::AU      !Number of kilometers per astronomical unit.
              REAL *8 ::CLIGHT  !Speed of light
              REAL *8 ::ASUN    !Solar radius
              REAL *8 ::RE      !Earth radius
              REAL *8 ::AM      !Moon radius
              CHARACTER*6 ::nam          !constant name be in use:TMP
              REAL *8 ::val     !constant value be in use:TMP
              
       !   Call constants
              nam='AU'
              CALL selconQ(nam,val)
              AU=val                   

              nam='CLIGHT'
              CALL selconQ(nam,val)     !km/s
              CLIGHT=val/AU             !AU/s

              nam='ASUN'
              CALL selconQ(nam,val)     !km
              ASUN=val/AU               !AU             

              nam='RE'
              CALL selconQ(nam,val)     !km
              RE=(val+65D0)/AU          !AU

              nam='AM' 
              CALL selconQ(nam,val)     !km
              AM=val/AU                 !AU
              
        !  Time conversion
              !Gregorian Calendar to Julian Calendar
              !Start_time
              DDay=(Second/3600+Minute/60+Hour)/24
              Day=Day+DDay
              CALL iau_CAL2JD ( Year, Month, Day, DJM0, DJM, J ) 
              
              !End_time
              DDay_1=(Second_1/3600+Minute_1/60+Hour_1)/24
              Day_1=Day_1+DDay_1              
              CALL iau_CAL2JD ( Year_1, Month_1, Day_1, DJM0_1, DJM_1, J )   
              
              
      !   Moon phase calculation
DO WHILE(DJM <= DJM_1)
    !Narrow the scope
    B=11 !Sun
    C=3  !Earth
    CALL PLEPH(DJM0+DJM,B,C,PO)
    E2S=PO(1:3)

    B=10 !Moon
    CALL PLEPH(DJM0+DJM,B,C,PO)
    E2M=PO(1:3)
    CALL iau_PDP(E2S,E2M,PRODUCT) 
    theta=ACOS(PRODUCT/(NORM2(E2S)*NORM2(E2M))) !(/rad)
    !WRITE(*,*)'theta',theta,'jd',DJM+DJM0
    IF(theta<THRESHOLD) THEN
            DJM=DJM+1
    ELSE
            DO
            !SUN 2 MOON 
            B=10 !Moon
            C=11 !Sun
            CALL PLEPH(DJM0+DJM,B,C,PO)           !TIME:Eclipse time
            S2M=PO(1:3)

            !Light time iteration(TO obtain the time when the Sun gives out light)
            dt0=0                                 !initial light time(/day)
            dt=NORM2(S2M)/CLIGHT/86400.0          !light time(/day)

            DO WHILE(ABS(dt0-dt)>light_time_diff)
            CALL PLEPH(DJM0+DJM-dt,B,C,PO)!TIME:the time when the Sun gives out light
            S2M=PO(1:3)
            dt0=dt
            dt=NORM2(S2M)/CLIGHT/86400.0      !update
            END DO

            t0=DJM0+DJM-dt                        !the time when the Sun gives out light
            B=11 !Sun
            C=3  !Earth
            CALL PLEPH(DJM0+DJM-dt,B,C,PO)        !TIME:the time when the Sun gives out light
            E2S=PO(1:3)
            !Light time iteration(To obtain the time when the light travel to earth)
            dt0=0!
            dt=NORM2(E2S)/CLIGHT/86400.0 

            DO WHILE(ABS(dt0-dt)>light_time_diff)
            CALL PLEPH(t0+dt,B,C,PO)               !TIME:the time when the light travel to earth
            E2S=PO(1:3)
            dt0=dt
            dt=NORM2(E2S)/CLIGHT/86400.0 
            END DO
            !Get the eclipse time	            
            B=10 
            C=11  
            CALL PLEPH(DJM0+DJM,B,C,PO)
            S2M=PO(1:3)
            E2O=RE*E2S/(RE-ASUN) 
            E2M=E2S+S2M    

            M2O=E2O-E2M 
            CALL iau_PDP(E2O,M2O,PRODUCT) 
            ANGLE_EOM=ACOS(PRODUCT/(NORM2(E2O)*NORM2(M2O))) 
            ANGLE_EO=ASIN(RE/NORM2(E2O)) 
            ANGLE_MO=ASIN(AM/NORM2(M2O)) 
            ERR=ANGLE_EOM-ANGLE_EO-ANGLE_MO 

            IF(NORM2(E2O)>NORM2(M2O).AND. ERR.LE.0D0) THEN
                    !Print the eclipse time
                    CALL iau_JD2CAL(DJM0,DJM,Year,Month,Day,DDay,J)    !
                    DDay=DDay-(32.184+37)*iteration_val
                    CALL sla_DD2TF(NDP,DDay,SIGN,IHMSF)!Convert an interval in days into hours, minutes, seconds
                    WRITE(*,*)'初亏时间：'
                    WRITE(*,*)Year,'年',Month,'月',Day,'日',IHMSF(1),':',IHMSF(2),':',IHMSF(3)+IHMSF(4)/1000.0  
                    DJM=DJM+1
                    EXIT   
            END IF
            DJM=DJM+iteration_val!iteration

            IF(DJM>DJM_1) EXIT

                  END DO
                  END IF
              END DO 
            
       END PROGRAM
