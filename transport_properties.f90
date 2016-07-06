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
    INTEGER :: iframe, ifile, jfile, jtype, jnum, inum, itype, iatom, itime, switch, tau
    REAL(4) :: drx,  dry,  drz
    REAL(4) :: drxp, dryp, drzp
    REAL(8) :: sumdrx, sumdry, sumdrz
    REAL(8) :: x_cm, x0_cm, y_cm, y0_cm, z_cm, z0_cm 

    ! mean-squared displacement in each direction (n_type,tsmax)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x0, y0, z0                    
    ! position at any given time origin (n_atoms(itype), t0max)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  x2t, y2t, z2t                
    ! number of steps between time origin
    INTEGER :: it0                                                                 
    ! counter for the sample branch --> "actual time"
    INTEGER :: ntel                                                                
    ! time origin
    INTEGER :: t0                                                                  
    ! corrected time origin; makes sure time origin number is not greater than t0max
    INTEGER :: tt0                                                                 
    ! difference in "actual time" and time origin
    INTEGER :: delta_t                                                             

    !!! UPDATE ACTUAL TIME FOR FRAME !!!
    ntel = ntel + 1            
 
    !!! START NEW TIME ORIGIN IF APPROPRIATE NUMBER OF TIMESTEPS HAVE PASSED !!!
    IF ((ntel .EQ. 1) .OR. (MOD(ntel, it0) .EQ. 0)) THEN
 
       !!! COUNT NUMBER OF TIME ORIGINS !!! 
       t0 = t0 + 1
 
       !!! STORE ONLY t0max TIME ORIGINS, COUNTER RESARTS IF EXCEEDED !!!
       tt0 = MOD(t0 - 1, t0max) + 1
 
       !!! RECORD STARTING TIME FOR EACH TIME ORIGIN !!!
       time0(tt0) = ntel
       DO iatom = 1, n_atoms_tot
          x0(iatom,tt0)  = rx(iatom)
          y0(iatom,tt0)  = ry(iatom)
          z0(iatom,tt0)  = rz(iatom)
       END DO
 
    END IF
 
    !!! STOP WHEN MAX NUMBER OF TIME ORIGINS IS REACHED !!!
    Tloop: DO tau = 1, MIN(t0,t0max)
    
       sumdrx = 0.0d0
       sumdry = 0.0d0
       sumdrz = 0.0d0
       x_cm  = 0.0d0
       y_cm  = 0.0d0
       z_cm  = 0.0d0
       x0_cm = 0.0d0
       y0_cm = 0.0d0
       z0_cm = 0.0d0
    
       !!! TIME STEP WITHIN TIME ORIGINS !!! 
       delta_t = ntel - time0(tau) + 1
       !!! STOP WHEN FINAL TIMESTEP IS REACHED!!!
       IF (delta_t .LT. tsmax) THEN
 
         ntime(delta_t) = ntime(delta_t) + 1
         itype = a_type(1)
         !!! LOOP OVER ALL ATOMS !!!
         Sloop: DO is = 1, nspecies
            IF (msd%t_skip(is) < 0) CYCLE
            Mloop: DO i = 1, nmolecules(is)
                im = locate(i, is)
                IF( .NOT. molecule_list(im, is)%live ) CYCLE
     
                Aloop: DO j = 1, 
                    
                !!! UPDATE MSDS !!!
                       
                drxp = (atom_list(ia, im, is)%rxp - x0(im, is, tau))
                dryp = (atom_list(ia, im, is)%ryp - y0(im, is, tau))
                drzp = (atom_list(ia, im, is)%rzp - z0(im, is, tau))

                CALL Minimum_Image_Separation(1, drxijp, dryijp, drzijp, drxij, dryij, drzij)
             
                x2t(itype,delta_t) = x2t(itype,delta_t) + (rx(iatom) - x0(iatom,tau))**2
                y2t(itype,delta_t) = y2t(itype,delta_t) + (ry(iatom) - y0(iatom,tau))**2
                z2t(itype,delta_t) = z2t(itype,delta_t) + (rz(iatom) - z0(iatom,tau))**2
                r2t(itype,delta_t) = (x2t(itype,delta_t) + y2t(itype,delta_t) + z2t(itype,delta_t))    
 
                x_cm  = x_cm + rx(iatom)
                y_cm  = y_cm + ry(iatom)
                z_cm  = z_cm + rz(iatom)
                x0_cm = x0_cm + x0(iatom,tau) 
                y0_cm = y0_cm + y0(iatom,tau) 
                z0_cm = z0_cm + z0(iatom,tau) 
       
                sumdrx = sumdrx + drx
                sumdry = sumdry + dry
                sumdrz = sumdrz + drz
     
                END DO Aloop
             END DO Mloop
         END DO Sloop
      END DO Tloop
     
  END SUBROUTINE Calculate_MSD

END MODULE Transport_Properties
