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

  SUBROUTINE Calculate_MSD(this_box)

    !*********************************************************************************
    !
    ! Calculate the Degree of ion Association
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: is, i, im, tau
    REAL(DP) :: drx,  dry,  drz
    REAL(DP) :: drxp, dryp, drzp

    ! difference in "actual time" and time origin
    INTEGER :: delta_t                                                             
    
    !!! UPDATE ACTUAL TIME FOR FRAME !!!
    msd%ntel = msd%ntel + 1            
 
    !!! START NEW TIME ORIGIN IF APPROPRIATE NUMBER OF TIMESTEPS HAVE PASSED !!!
    !IF ((msd%ntel .EQ. 1) .OR. (MOD(msd%ntel, msd%t_origin) .EQ. 0)) THEN
    IF (mod(msd%ntel-1, msd%t_origin) .EQ. 0) THEN
 
        !!! COUNT NUMBER OF TIME ORIGINS !!! 
        msd%t0 = msd%t0 + 1
  
        !!! STORE ONLY msd%t0max TIME ORIGINS, COUNTER RESARTS IF EXCEEDED !!!
        msd%tt0 = MOD(msd%t0 - 1, msd%t0max) + 1
  
        !!! RECORD STARTING TIME FOR EACH TIME ORIGIN !!!
        msd%time0(msd%tt0) = msd%ntel
        DO is = 1, nspecies
            IF (.not. msd%species(is)) CYCLE

            DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE

                msd%x0(im, is, msd%tt0) = molecule_list(im, is)%xcom
                msd%y0(im, is, msd%tt0) = molecule_list(im, is)%ycom
                msd%z0(im, is, msd%tt0) = molecule_list(im, is)%zcom
            END DO
        END DO
 
    END IF
 
    !!! STOP WHEN MAX NUMBER OF TIME ORIGINS IS REACHED !!!
    Tloop: DO tau = 1, MIN(msd%t0, msd%t0max)
    
        !!! TIME STEP WITHIN TIME ORIGINS !!! 
        delta_t = msd%ntel - msd%time0(tau) + 1
        !!! STOP WHEN FINAL TIMESTEP IS REACHED!!!
        !IF (delta_t .GT. n_mcsteps) CYCLE
  
        msd%ntime(delta_t) = msd%ntime(delta_t) + 1
        !!! LOOP OVER ALL ATOMS !!!
        Sloop: DO is = 1, nspecies
            IF (.not. msd%species(is)) CYCLE

            Mloop: DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
                    
                !!! UPDATE MSDS !!!
                drx = molecule_list(im, is)%xcom - msd%x0(im, is, tau)
                dry = molecule_list(im, is)%ycom - msd%y0(im, is, tau)
                drz = molecule_list(im, is)%zcom - msd%z0(im, is, tau)
 
                ! this needs to be from unwrapped trajcetories
                !CALL Minimum_Image_Separation(1, drxp, dryp, drzp, drx, dry, drz)
             
                msd%x2t(is,delta_t) = msd%x2t(is,delta_t) + (drx * drx)
                msd%y2t(is,delta_t) = msd%y2t(is,delta_t) + (dry * dry)
                msd%z2t(is,delta_t) = msd%z2t(is,delta_t) + (drz * drz)
                msd%r2t(is,delta_t) = msd%r2t(is,delta_t) + ((drx * drx) + (dry * dry) + (drz * drz))
 
            END DO Mloop
        END DO Sloop
    END DO Tloop
     
  END SUBROUTINE Calculate_MSD

END MODULE Transport_Properties
