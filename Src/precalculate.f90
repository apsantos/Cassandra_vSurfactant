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

SUBROUTINE precalculate
!*******************************************************************************  
    ! Revision history:
    !
    ! 12/10/13  : Beta version 
!********************************************************************************    
    USE Run_variables
    USE File_Names
    USE Energy_Routines
    IMPLICIT NONE
    
    
    INTEGER :: ibox, itype, jtype, is, ia, ja
    ! Determine the direct sum maximum cutoff distances squared, beyond 
    ! which we ingore pairwise interactions. 
    ! Will need to alter when neighborlist added.
    
    REAL(DP) :: roffsq_ronsq

    roffsq_ronsq = 0.0_DP

    DO ibox = 1,nbr_boxes
       IF (int_vdw_style(ibox) /= vdw_none) THEN
          IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
             ! Compute the square for use in pair routines
             ron_charmmsq(ibox) = ron_charmm(ibox) * ron_charmm(ibox)
             roff_charmmsq(ibox) = roff_charmm(ibox) * roff_charmm(ibox)

          ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN

             ron_switch_sq(ibox) = ron_switch(ibox) * ron_switch(ibox)

             roff_switch_sq(ibox) = roff_switch(ibox) * roff_switch(ibox)

             roffsq_ronsq = roff_switch_sq(ibox) - ron_switch_sq(ibox)

             switch_factor1(ibox) = roffsq_ronsq**3.0_DP
             switch_factor1(ibox) = 1.0_DP / switch_factor1(ibox)

             switch_factor2(ibox) = roff_switch_sq(ibox) - (3.0_DP * ron_switch_sq(ibox))

             rcut_vdw(ibox) = roff_switch(ibox)

          
          ELSE
             ! Compute the square for use in pair routines
             DO itype = 1, nbr_atomtypes
                DO jtype = 1, nbr_atomtypes
                    rcut_vdwsq_mix(itype, jtype) = rcut_vdw_mix(itype, jtype) * rcut_vdw_mix(itype, jtype)
                END DO
             END DO
             DO is = 1, nspecies
                DO ia = 1, natoms(is)
                    DO ja = 1, natoms(is)
                       rcut_in_vdwsq_mix(ia, ja, is) = rcut_in_vdw_mix(ia, ja, is) * rcut_in_vdw_mix(ia, ja, is)
                    END DO
                END DO
             END DO

             rcut_vdwsq(ibox) = rcut_vdw(ibox) * rcut_vdw(ibox)
             rcut_vdw3(ibox) = rcut_vdwsq(ibox) * rcut_vdw(ibox)
             rcut_vdw6(ibox) = rcut_vdw3(ibox) * rcut_vdw3(ibox)
          ENDIF

       ENDIF
    
       IF (int_charge_style(ibox) /= charge_none) THEN
          rcut_coulsq(ibox) = rcut_coul(ibox) * rcut_coul(ibox)
       ENDIF
    
       IF ( (int_vdw_style(ibox) /= vdw_none) .AND. (int_charge_style(ibox) /= charge_none)) THEN
          rcut_max(ibox) = MAX(rcut_vdw(ibox),rcut_coul(ibox))
       ELSE IF ( int_vdw_style(ibox) /= vdw_none) THEN
          rcut_max(ibox) = rcut_vdw(ibox)
       ELSE
          rcut_max(ibox) = rcut_coul(ibox)
       END IF

    END DO

    rcut_lowsq = rcut_low * rcut_low

    ! initialize the pair_nrg array
    
    IF (l_pair_nrg) THEN
       pair_nrg_vdw(:,:) = 0.0_DP
       pair_nrg_qq(:,:) = 0.0_DP
    END IF

    ! Add other stuff like molecule mass, LRC, and anything else needed.
    ! ALLOCATE memory for the Ewald stuff

    ALLOCATE(energy(nbr_boxes),virial(nbr_boxes))
    energy(:)%inter_vdw = 0.0_DP
    energy(:)%lrc = 0.0_DP
    energy(:)%inter_q = 0.0_DP
    energy(:)%intra_vdw = 0.0_DP
    energy(:)%intra_q = 0.0_DP
    energy(:)%intra = 0.0_DP
    energy(:)%ewald_reciprocal = 0.0_DP
    energy(:)%ewald_self = 0.0_DP
    energy(:)%total = 0.0_DP
    energy(:)%bond = 0.0_DP
    energy(:)%angle = 0.0_DP
    energy(:)%dihedral = 0.0_DP
    energy(:)%improper = 0.0_DP
    energy(:)%erf_self = 0.0_DP
    energy(:)%ewald_self_calc = .false.

    virial(:)%inter_vdw = 0.0_DP
    virial(:)%lrc = 0.0_DP
    virial(:)%inter_q = 0.0_DP
    virial(:)%intra_vdw = 0.0_DP
    virial(:)%intra_q = 0.0_DP
    virial(:)%intra = 0.0_DP
    virial(:)%ewald_reciprocal = 0.0_DP
    virial(:)%ewald_self = 0.0_DP
    virial(:)%total = 0.0_DP
    virial(:)%bond = 0.0_DP
    virial(:)%angle = 0.0_DP
    virial(:)%dihedral = 0.0_DP
    virial(:)%improper = 0.0_DP
    virial(:)%erf_self = 0.0_DP
    virial(:)%ewald_self_calc = .false.

  END SUBROUTINE precalculate
