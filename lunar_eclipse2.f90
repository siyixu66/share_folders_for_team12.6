SUBROUTINE compute_angle(DJM0,DJM,ASUN,RE,AM,CLIGHT,ANGLE_EOM,ANGLE_EO,ANGLE_MO,E2O,M2O)
        REAL *8 :: DJM0,DJM                       !输入儒略日
        REAL *8,DIMENSION(6) :: PO                !输出坐标
        INTEGER :: B                              !目标天体编号
        INTEGER :: C                              !中心天体编号
        REAL *8::dt0,dt  
        REAL *8,DIMENSION(3)::E2S                 !地日距离
        REAL *8,DIMENSION(3)::E2O                 !地心与影锥距离
        REAL *8,DIMENSION(3)::E2M               !地月距离
        REAL *8,DIMENSION(3)::M2O                 !月心与影锥的距离
        REAL *8,DIMENSION(3)::S2M                 !日月矢量
        REAL *8::ANGLE_EOM                        !地锥月角度
        REAL *8::ANGLE_EO                         !地锥角
        REAL *8::ANGLE_MO                         !月锥角
        REAL *8::PRODUCT                          !内积

        REAL*8::light_time_diff=0.001/86400.0        !两次迭代的光行时差，单位为天

        
        REAL *8::ASUN                             !太阳半径
        REAL *8::RE                               !地球赤道半径
        REAL *8::AM                               !月球半径
        REAL *8::CLIGHT                           !光速
      
        B=11 !sun
        C=3  !earth
        CALL PLEPH(DJM0+DJM,B,C,PO)!太阳相对于地球的位置
        E2S=PO(1:3)
        !光行时迭代
        dt=0!光行时，单位为天
        dt0=NORM2(E2S)/CLIGHT/86400.0 !

        DO WHILE(ABS(dt0-dt)>light_time_diff)
            dt0=dt
            CALL PLEPH(DJM0+DJM-dt0,B,C,PO)!太阳相对于地球的位置,单位为AU
            E2S=PO(1:3)
            dt=NORM2(E2S)/CLIGHT/86400.0 !光行时，单位为天
        END DO

        E2O=RE*E2S/(RE-ASUN) !地心与影锥距离
        B=10 !moon
        CALL PLEPH(DJM0+DJM,B,C,PO)!月球相对于地球的位置
        E2M=PO(1:3)
        M2O=E2O-E2M !月心到影锥
        CALL iau_PDP(E2O,M2O,PRODUCT) !内积
        ANGLE_EOM=ACOS(PRODUCT/(NORM2(E2O)*NORM2(M2O))) !地锥月角度
        ANGLE_EO=ASIN(RE/NORM2(E2O)) !地锥角
        ANGLE_MO=ASIN(AM/NORM2(M2O)) !月锥角
END SUBROUTINE compute_angle

SUBROUTINE print_time(DJM0,DJM,iteration_val,name)
        REAL *8 :: DJM0,DJM                       !输入儒略日
        INTEGER :: YEAR,MONTH,DAY,J               !输出年月日
        REAL *8 :: DDAY                            !日小数部分
        INTEGER :: HOUR                           !时
        INTEGER :: MIN                            !分
        REAL *8 :: SECOND                          !秒
        INTEGER :: NDP=3                          !小数秒位数
        INTEGER ,DIMENSION(4) :: IHMSF            !时分秒小数秒数组
        CHARACTER :: SIGN                         !正负号
        CHARACTER*12 :: name                       !常数名字
        REAL *8::iteration_val
        CALL iau_JD2CAL(DJM0,DJM,YEAR,MONTH,DAY,DDAY,J)
        DDAY=DDAY-(32.184+37)*iteration_val
        CALL sla_DD2TF(NDP,DDAY,SIGN,IHMSF)!小数天数转为时分秒
        WRITE(*,*)name
        WRITE(*,*)YEAR,'年',MONTH,'月',DAY,'日',IHMSF(1),':',IHMSF(2),':',IHMSF(3)+IHMSF(4)/1000.0
END SUBROUTINE print_time

PROGRAM MOON
          implicit none
          REAL *8,DIMENSION(6) :: PO                !输出坐标
          INTEGER :: B                              !目标天体编号
          INTEGER :: C                              !中心天体编号    
          INTEGER :: YEAR=2025                     !年
          INTEGER :: MONTH=1                        !月
          INTEGER :: DAY=1                          !日

          INTEGER::YEAR1=2026                     !结束的年
          INTEGER::MONTH1=1                         !结束的月
          INTEGER::DAY1=1                           !结束的日

          REAL*8 :: DDAY                            !日小数部分
          INTEGER :: HOUR                           !时
          INTEGER :: MIN                            !分
          REAL*8 :: SECOND                          !秒
          REAL *8 :: DJM0                           !儒略日1
          REAL *8 :: DJM                            !儒略日2

          REAL*8::DJM0_1                            !结束的儒略日1
          REAL*8::DJM_1                             !结束的儒略日2

          REAL*8 TEMP1
          REAL*8 TEMP2 
          INTEGER :: J                              !状态
          INTEGER :: NDP=3                          !小数秒位数
          INTEGER ,DIMENSION(4) :: IHMSF            !时分秒小数秒数组
          CHARACTER :: SIGN                         !正负号
          CHARACTER*6 :: nams                       !常数名字
          REAL *8::vals                             !常数值
          REAL *8::iteration_val=1.0/(24*3600)      !循环迭代量，1s为多少天
          REAL*8::light_time_diff=0.001/86400.0        !两次迭代的光行时差，单位为天

          REAL *8::ERR                              !角度差，判断是否发生食既
          REAL *8::ASUN                             !太阳半径
          REAL *8::RE                               !地球赤道半径
          REAL *8::AM                               !月球半径
          REAL *8::CLIGHT                           !光速
          REAL *8::AU                               !天文单位          
          REAL *8,DIMENSION(3)::E2S                 !地日距离
          
          REAL *8,DIMENSION(3)::E2O                 !地心与影锥距离
          REAL *8,DIMENSION(3)::E2M               !地月距离

          REAL*8::theta                             !地月夹角
          REAL *8,DIMENSION(3)::M2O                 !月心与影锥的距离
          REAL *8,DIMENSION(3)::S2M                 !日月矢量
          REAL *8::ANGLE_EOM                        !地锥月角度
          REAL*8::ANGLE_EOM1
          REAL *8::ANGLE_EO                         !地锥角
          REAL *8::ANGLE_MO                          !月锥角
          REAL*8::delta                        !食甚的时候月球距离阴影面中心的距离
          REAL*8::MAG                           !食分
          REAL*8::SHADOW_RADIUS                 !阴影半径
          REAL *8::PRODUCT                          !内积

          
          REAL *8, PARAMETER :: PI = 3.14159265358979323846
          REAL *8, PARAMETER :: THRESHOLD = 175.0 * PI / 180.0
        
          REAL *8::dt0,dt                          !光行时
!**************读出常量***************
          nams='AU'         !天文单位
          CALL selconQ(nams,vals)
          AU=VALS          
          nams='ASUN'       !太阳半径
          CALL selconQ(nams,vals)
          ASUN=VALS/AU          
          nams='RE'         !地球半径(atmosphere)
          CALL selconQ(nams,vals)
          RE=(VALS+65D0)/AU        
          nams='AM'         !月球半径
          CALL selconQ(nams,vals)
          AM=VALS/AU          
          nams='CLIGHT'     !光速
          CALL selconQ(nams,vals)
          CLIGHT=VALS/AU !光速，单位为AU/s
          
          CALL iau_CAL2JD ( YEAR, MONTH, DAY, DJM0, DJM, J )    !格里历转为儒略日
          CALL iau_CAL2JD ( YEAR1, MONTH1, DAY1, DJM0_1, DJM_1, J )    !格里历转为儒略日
          
!往后递推，只看日地月夹角大于175度
DO WHILE(DJM <= DJM_1)
    B=11 !sun
    C=3  !earth
    CALL PLEPH(DJM0+DJM,B,C,PO)!太阳相对于地球的位置
    E2S=PO(1:3)
    B=10 !moon
    CALL PLEPH(DJM0+DJM,B,C,PO)!月球相对于地球的位置
    E2M=PO(1:3)
    CALL iau_PDP(E2S,E2M,PRODUCT) !内积
    theta=ACOS(PRODUCT/(NORM2(E2S)*NORM2(E2M))) !地月夹角,rad
    
    IF(theta<THRESHOLD) THEN
        DJM=DJM+1   
    ELSE
        !按照1s的迭代量，往后寻找月食初亏时刻   
    DO 
        CALL compute_angle(DJM0,DJM,ASUN,RE,AM,CLIGHT,ANGLE_EOM,ANGLE_EO,ANGLE_MO,E2O,M2O)
        ERR=ANGLE_EOM-ANGLE_EO-ANGLE_MO !角度差
        IF(NORM2(E2O)>NORM2(M2O).AND. ERR.LE.0D0) THEN
        !月食初亏时刻输出
            CALL print_time(DJM0,DJM,iteration_val,'初亏时间：')
            !寻找食甚时刻，即地锥月角度最小     
            CALL compute_angle(DJM0,DJM+iteration_val,ASUN,RE,AM,CLIGHT,ANGLE_EOM1,ANGLE_EO,ANGLE_MO,E2O,M2O)
            DO WHILE(ANGLE_EOM1<ANGLE_EOM)
                ANGLE_EOM=ANGLE_EOM1
                DJM=DJM+iteration_val
                CALL compute_angle(DJM0,DJM+iteration_val,ASUN,RE,AM,CLIGHT,ANGLE_EOM1,ANGLE_EO,ANGLE_MO,E2O,M2O)
            END DO
            CALL print_time(DJM0,DJM,iteration_val,'食甚时间：')
            !计算食分
            delta=NORM2(M2O)*SIN(ANGLE_EOM)
            SHADOW_RADIUS=RE*NORM2(M2O)/NORM2(E2O)
            MAG=((AM-delta)+SHADOW_RADIUS)/(2*AM)
            WRITE(*,*)'食分：',MAG
            IF(MAG>=1)THEN
                !食既，仅月全食有。月球内切阴影
                TEMP1=DJM
                TEMP2=ANGLE_EOM
                DO WHILE(ANGLE_EOM<ANGLE_EO-ANGLE_MO)
                    CALL compute_angle(DJM0,DJM,ASUN,RE,AM,CLIGHT,ANGLE_EOM,ANGLE_EO,ANGLE_MO,E2O,M2O)
                    DJM=DJM-iteration_val
                END DO
                CALL print_time(DJM0,DJM,iteration_val,'食既时间：')
                DJM=TEMP1
                ANGLE_EOM=TEMP2
                !计算生光时刻，只有月全食才有
                DO WHILE(ANGLE_EOM<ANGLE_EO-ANGLE_MO)
                    CALL compute_angle(DJM0,DJM,ASUN,RE,AM,CLIGHT,ANGLE_EOM,ANGLE_EO,ANGLE_MO,E2O,M2O)
                    DJM=DJM+iteration_val
                END DO
                CALL print_time(DJM0,DJM,iteration_val,'生光时间：')
            END IF
            
        !寻找月食复圆时刻，即再次外切时刻
            DO WHILE(ANGLE_EOM<ANGLE_EO+ANGLE_MO)
                CALL compute_angle(DJM0,DJM,ASUN,RE,AM,CLIGHT,ANGLE_EOM,ANGLE_EO,ANGLE_MO,E2O,M2O)
                DJM=DJM+iteration_val
            END DO
            CALL print_time(DJM0,DJM,iteration_val,'复圆时间：')

            DJM=DJM+1
            EXIT   
        END IF
        DJM=DJM+iteration_val!儒略日+1s

        IF(DJM>DJM_1) EXIT

    END DO
    END IF  
END DO 
END PROGRAM MOON
