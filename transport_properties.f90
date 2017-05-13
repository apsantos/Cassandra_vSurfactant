!********************************************************************************
!   Cassandra - An open source atomistic Monte Carlo software package
!   developed at the University of Notre Dame.
!   http://cassandra.nd.edu
!   Prof. Edward Maginn <ed@nd.edu>
!   Copyright (2013) University of Notre Dame du Lac
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!********************************************************************************

MODULE Transport_Properties

  !******************************************************************************
  ! 
  ! - Andrew P. Santos
  !
  !*******************************************************************************

  USE Run_Variables
  USE IO_Utilities
  USE Cluster_Routines

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Calculate_MSD()

    !*********************************************************************************
    !
    ! Calculate the Mean-squared Displacement
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER :: is, i, im, tau
    REAL(DP) :: drx,  dry,  drz

    ! difference in "actual time" and time origin
    INTEGER :: delta_t                                                             
    
    !!! UPDATE ACTUAL TIME FOR FRAME !!!
    trans%ntel = trans%ntel + 1            
 
    !!! START NEW TIME ORIGIN IF APPROPRIATE NUMBER OF TIMESTEPS HAVE PASSED !!!
    !IF ((trans%ntel .EQ. 1) .OR. (MOD(trans%ntel, trans%t_origin) .EQ. 0)) THEN
    IF (mod(trans%ntel-1, trans%t_origin) .EQ. 0) THEN
 
        !!! COUNT NUMBER OF TIME ORIGINS !!! 
        trans%t0 = trans%t0 + 1
  
        !!! STORE ONLY trans%t0max TIME ORIGINS, COUNTER RESARTS IF EXCEEDED !!!
        trans%tt0 = MOD(trans%t0 - 1, trans%t0max) + 1
  
        !!! RECORD STARTING TIME FOR EACH TIME ORIGIN !!!
        trans%time0(trans%tt0) = trans%ntel
        DO is = 1, nspecies
            IF (.not. trans%msd_species(is)) CYCLE

            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE

                trans%rx0(im, is, trans%tt0) = molecule_list(im, is)%xcom
                trans%ry0(im, is, trans%tt0) = molecule_list(im, is)%ycom
                trans%rz0(im, is, trans%tt0) = molecule_list(im, is)%zcom
            END DO
        END DO
 
    END IF
 
    !!! STOP WHEN MAX NUMBER OF TIME ORIGINS IS REACHED !!!
    Tloop: DO tau = 1, MIN(trans%t0, trans%t0max)
    
        !!! TIME STEP WITHIN TIME ORIGINS !!! 
        delta_t = trans%ntel - trans%time0(tau) + 1
        !!! STOP WHEN FINAL TIMESTEP IS REACHED!!!
        !IF (delta_t .GT. n_mcsteps) CYCLE
  
        trans%ntime(delta_t) = trans%ntime(delta_t) + 1
        !!! LOOP OVER ALL ATOMS !!!
        Sloop: DO is = 1, nspecies
            IF (.not. trans%msd_species(is)) CYCLE

            Mloop: DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                    
                !!! UPDATE MSDS !!!
                drx = molecule_list(im, is)%xcom - trans%rx0(im, is, tau)
                dry = molecule_list(im, is)%ycom - trans%ry0(im, is, tau)
                drz = molecule_list(im, is)%zcom - trans%rz0(im, is, tau)
 
                ! this needs to be from unwrapped trajcetories
                !CALL Minimum_Image_Separation(1, drxp, dryp, drzp, drx, dry, drz)
             
                trans%x_msd(is,delta_t) = trans%x_msd(is,delta_t) + (drx * drx)
                trans%y_msd(is,delta_t) = trans%y_msd(is,delta_t) + (dry * dry)
                trans%z_msd(is,delta_t) = trans%z_msd(is,delta_t) + (drz * drz)
                !trans%msd(is,delta_t) = trans%msd(is,delta_t) + ((drx * drx) + (dry * dry) + (drz * drz))
 
            END DO Mloop
        END DO Sloop
    END DO Tloop
     
  END SUBROUTINE Calculate_MSD

  SUBROUTINE Calculate_VACF()

    !*********************************************************************************
    !
    ! Calculate the Velocity Autocorrelation function
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER :: is, i, im, tau
    !REAL(DP) :: drx,  dry,  drz

    ! difference in "actual time" and time origin
    INTEGER :: delta_t                                                             
    
    !!! UPDATE ACTUAL TIME FOR FRAME !!!
    trans%ntel = trans%ntel + 1            
 
    !!! START NEW TIME ORIGIN IF APPROPRIATE NUMBER OF TIMESTEPS HAVE PASSED !!!
    !IF ((trans%ntel .EQ. 1) .OR. (MOD(trans%ntel, trans%t_origin) .EQ. 0)) THEN
    IF (mod(trans%ntel-1, trans%t_origin) .EQ. 0) THEN
 
        !!! COUNT NUMBER OF TIME ORIGINS !!! 
        trans%t0 = trans%t0 + 1
  
        !!! STORE ONLY trans%t0max TIME ORIGINS, COUNTER RESARTS IF EXCEEDED !!!
        trans%tt0 = MOD(trans%t0 - 1, trans%t0max) + 1
  
        !!! RECORD STARTING TIME FOR EACH TIME ORIGIN !!!
        trans%time0(trans%tt0) = trans%ntel
        DO is = 1, nspecies
            IF (.not. trans%vacf_species(is)) CYCLE

            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE

                trans%vx0(im, is, trans%tt0) = molecule_list(im, is)%vxcom
                trans%vy0(im, is, trans%tt0) = molecule_list(im, is)%vycom
                trans%vz0(im, is, trans%tt0) = molecule_list(im, is)%vzcom
            END DO
        END DO
 
    END IF
 
    !!! STOP WHEN MAX NUMBER OF TIME ORIGINS IS REACHED !!!
    Tloop: DO tau = 1, MIN(trans%t0, trans%t0max)
    
        !!! TIME STEP WITHIN TIME ORIGINS !!! 
        delta_t = trans%ntel - trans%time0(tau) + 1
        !!! STOP WHEN FINAL TIMESTEP IS REACHED!!!
        !IF (delta_t .GT. n_mcsteps) CYCLE
  
        trans%ntime(delta_t) = trans%ntime(delta_t) + 1
        !!! LOOP OVER ALL ATOMS !!!
        Sloop: DO is = 1, nspecies
            IF (.not. trans%vacf_species(is)) CYCLE

            Mloop: DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                    
                !!! UPDATE MSDS !!!
                ! this needs to be from unwrapped trajcetories
                !CALL Minimum_Image_Separation(1, drxp, dryp, drzp, drx, dry, drz)
             
                trans%x_vacf(is,delta_t) = trans%x_vacf(is,delta_t) + (molecule_list(im, is)%vxcom * trans%vx0(im, is, tau))
                trans%y_vacf(is,delta_t) = trans%y_vacf(is,delta_t) + (molecule_list(im, is)%vycom * trans%vy0(im, is, tau))
                trans%z_vacf(is,delta_t) = trans%z_vacf(is,delta_t) + (molecule_list(im, is)%vzcom * trans%vz0(im, is, tau))
                trans%vacf(is,delta_t) = trans%vacf(is,delta_t) + &
                                         (molecule_list(im, is)%vxcom * trans%vx0(im, is, tau)) + &
                                         (molecule_list(im, is)%vycom * trans%vy0(im, is, tau)) + &
                                         (molecule_list(im, is)%vzcom * trans%vz0(im, is, tau)) 
 
            END DO Mloop
        END DO Sloop
    END DO Tloop
     
  END SUBROUTINE Calculate_VACF

END MODULE Transport_Properties
