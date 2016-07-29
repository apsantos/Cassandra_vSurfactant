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

MODULE End_To_End

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

  SUBROUTINE Calculate_End_To_End_Distance(this_box)

    !*********************************************************************************
    !
    ! Calculate the End_To_End Distance of species as a function of cluster size
    !
    ! 2/19/15  : Andrew P. Santos
    !*********************************************************************************

    INTEGER, INTENT(IN) :: this_box
    INTEGER :: i, j, is, im
    INTEGER :: clus_size
    REAL(DP) :: rxij, ryij, rzij, rijsq, rxijp, ryijp, rzijp

    DO is = 1, nspecies
        IF (.not. end2end%species(is) ) CYCLE

        DO i = 1, nmolecules(is)
            im = locate(i, is)
            IF( .NOT. molecule_list(im, is)%live ) CYCLE

            ! Get the positions of the COM of the two molecule species
            rxijp = atom_list(1, im, is)%rxp - atom_list(natoms(is), im, is)%rxp
            ryijp = atom_list(1, im, is)%ryp - atom_list(natoms(is), im, is)%ryp
            rzijp = atom_list(1, im, is)%rzp - atom_list(natoms(is), im, is)%rzp
            
            ! Now get the minimum image separation 
            CALL Minimum_Image_Separation(1, rxijp, ryijp, rzijp, rxij, ryij, rzij)
    
            rijsq = (rxij*rxij) + (ryij*ryij) + (rzij*rzij)
            
            clus_size = cluster%N( cluster%clabel(im, is) )
            end2end%distance(clus_size, is) = end2end%distance(clus_size, is) + ( SQRT(rijsq) / DBLE(clus_size) )

        END DO
    END DO
                    
  END SUBROUTINE Calculate_End_To_End_Distance

END MODULE End_To_End
