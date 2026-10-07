module yog_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Yoval-O'Gorman deep convection scheme.
!   based off of the Zhang-McFarlane deep convection scheme
!
! Author: J.W. Atkinson
! January 2024
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use cam_abortutils,   only: endrun
   use physconst,        only: pi
   use spmd_utils,       only: masterproc
   use perf_mod
   use cam_logfile,  only: iulog
   use constituents, only: cnst_add

   implicit none
   private
   save

   ! Public methods

   public ::&
      yog_register,           & 
      yog_readnl,             &! read namelist for YOG scheme/module
      yog_init,               &! initialize YOG scheme/module
      yog_final,              &! finalize YOG scheme/module
      yog_tend                 ! return tendencies

   character(len=136) :: yog_nn_weights    ! location of weights for the YOG NN, set in namelist
   character(len=136) :: SAM_sounding  ! location of SAM sounding profile for the YOG NN, set in namelist

   integer :: prec_dp_idx, snow_dp_idx
   integer :: dp_flxprc_idx, dp_flxsnw_idx, dp_cldice_idx, dp_cldliq_idx

!=========================================================================================
contains
!=========================================================================================

! There is currently no need for a yog_register as nothing to add to the physics buffer
subroutine yog_register

 !----------------------------------------
 ! Purpose: register fields with the physics buffer
 !----------------------------------------

   use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4

   implicit none

   integer idx

   ! XL XL XL
   !call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx)
   !call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx)
   !call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)
   !call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)   
   
end subroutine yog_register

!=========================================================================================

subroutine yog_readnl(nlfile)

   use spmd_utils,      only: mpi_character, masterprocid, mpicom
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'yog_readnl'

   namelist /yog_params_nl/ yog_nn_weights, SAM_sounding
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'yog_params_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, yog_params_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

   end if

   ! Broadcast namelist variables
   call mpi_bcast(yog_nn_weights,   len(yog_nn_weights),   mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("yog_readnl: FATAL: mpi_bcast: yog_nn_weights")
   call mpi_bcast(SAM_sounding, len(SAM_sounding), mpi_character, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("yog_readnl: FATAL: mpi_bcast: SAM_sounding")

end subroutine yog_readnl

!=========================================================================================

subroutine yog_init()

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld, add_default, horiz_only
  use nn_interface_CAM,   only: nn_convection_flux_CAM_init

  use spmd_utils,     only: masterproc
  use phys_control,   only: phys_getopts

  integer istat
  logical :: history_budget ! output tendencies and state variables for
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

  ! Register fields with the output buffer
  call addfld ('YOGDT  ',   (/ 'lev' /),  'A', 'K/s','T tendency - Yuval-OGorman moist convection')
  call addfld ('YOGDQ  ',   (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Yuval-OGorman moist convection')
  call addfld ('YOGDICE',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud ice tendency - Yuval-OGorman convection')
  call addfld ('YOGDLIQ',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud liq tendency - Yuval-OGorman convection')
  call addfld ('PREC_YOG',  horiz_only ,  'A', 'm/s','Surface preciptation - Yuval-OGorman convection')
  call addfld ('SNOW_YOG',  horiz_only ,  'A', 'm/s','Surface snow - Yuval-OGorman convection')  
  call addfld ('YOGDNUMLIQ', (/ 'lev' /),  'A', 'N/s','Cloud liq number conc. tendency - Yuval-OGorman convection')
  call addfld ('YOGDNUMICE', (/ 'lev' /),  'A', 'N/s','Cloud ice number conc. tendency - Yuval-OGorman convection')
  
  call addfld ('YOG_DQ_ADV  ',   (/ 'lev' /),  'A', '(kg/kg)/kg/s','Q subg. vert. adv. flux conv. - YOG ')
  call addfld ('YOG_DQ_AUTO ',   (/ 'lev' /),  'A', '(kg/kg)/kg/s','Q cloud microphysics tend. - YOG')
  call addfld ('YOG_DQ_SED  ',   (/ 'lev' /),  'A', '(kg/kg)/kg/s','Q cloud ice sediment. flux conv. - YOG')
  call addfld ('YOG_DT_ADV  ',   (/ 'lev' /),  'A', 'J/kg/s','H subg. vert. adv. flux conv. - YOG')
  call addfld ('YOG_DT_AUTO ',   (/ 'lev' /),  'A', 'J/kg/s','H cloud microphysics tend. - YOG')
  call addfld ('YOG_DT_SED  ',   (/ 'lev' /),  'A', 'J/kg/s','H cloud ice sediment. flux conv. - YOG')
  call addfld ('YOG_DT_RAD  ',   (/ 'lev' /),  'A', 'J/kg/s','H radiative & phase change tend. - YOG')

  call addfld ('YOG_DT_ADVFLX  ',   (/ 'lev' /),  'A', 'J/m2/s','H subg. vert. adv. flux - YOG')
  
!  call addfld ('DP_CLDLIQ', (/ 'lev' /),  'A', 'kg/kg', 'deep gbm cloud liquid water')
!  call addfld ('DP_CLDICE', (/ 'lev' /),  'A', 'kg/kg', 'deep gbm cloud ice water')
!  call addfld ('DP_FLXPRC', (/ 'lev' /),  'A', 'kg/m2/s', 'deep interface gbm convective precipitation flux')
!  call addfld ('DP_FLXSNW', (/ 'lev' /),  'A', 'kg/m2/s', 'deep interface gbm convective snow flux')    
  
  if (masterproc) then
     write(iulog,*)'YOG output fields added to buffer'
  end if

  call phys_getopts( history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num)

  if ( history_budget ) then
     call add_default('YOGDT  ' , history_budget_histfile_num, ' ')
     call add_default('YOGDQ  ' , history_budget_histfile_num, ' ')
     call add_default('YOGDICE' , history_budget_histfile_num, ' ')
     call add_default('YOGDLIQ' , history_budget_histfile_num, ' ')
     call add_default('PREC_YOG', history_budget_histfile_num, ' ')
     call add_default('SNOW_YOG', history_budget_histfile_num, ' ')     
     call add_default('YOGDNUMLIQ', history_budget_histfile_num, ' ')
     call add_default('YOGDNUMICE', history_budget_histfile_num, ' ')     
     
     call add_default('YOG_DQ_ADV ' , history_budget_histfile_num, ' ')
     call add_default('YOG_DQ_AUTO' , history_budget_histfile_num, ' ')
     call add_default('YOG_DQ_SED ' , history_budget_histfile_num, ' ')
     call add_default('YOG_DT_ADV ' , history_budget_histfile_num, ' ')
     call add_default('YOG_DT_AUTO' , history_budget_histfile_num, ' ')
     call add_default('YOG_DT_SED ' , history_budget_histfile_num, ' ')
     call add_default('YOG_DT_RAD ' , history_budget_histfile_num, ' ')           

     call add_default('YOG_DT_ADVFLX ' , history_budget_histfile_num, ' ')     

     ! XL XL XL     
     !call add_default('DP_CLDLIQ', history_budget_histfile_num, ' ')
     !call add_default('DP_CLDICE', history_budget_histfile_num, ' ')
     
  end if

  call nn_convection_flux_CAM_init(yog_nn_weights, SAM_sounding)
  if (masterproc) then
     write(iulog,*)'yog_nn_weights at: ', yog_nn_weights
     write(iulog,*)'SAM_sounding at: ', SAM_sounding
     write(iulog,*)'YOG scheme initialised'
  endif

end subroutine yog_init

!=========================================================================================

subroutine yog_final()

!----------------------------------------
! Purpose:  finalization of the YOG scheme
!----------------------------------------

  use nn_interface_CAM,   only: nn_convection_flux_CAM_finalize

  call nn_convection_flux_CAM_finalize()

end subroutine yog_final

!=========================================================================================

subroutine yog_tend(ztodt, state, ptend, pbuf, landfrac, ideep)

!----------------------------------------
! Purpose:  tendency calculation for YOG scheme
!----------------------------------------

   use cam_history,    only: outfld
   use physics_types,  only: physics_state, physics_ptend
   use physics_types,  only: physics_ptend_init
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   use constituents,   only: pcnst, cnst_get_ind
   
   use nn_interface_CAM,   only: nn_convection_flux_CAM

   ! Arguments

   type(physics_state), intent(in)          :: state          ! Physics state variables
   type(physics_ptend), intent(out)         :: ptend          ! individual parameterization tendencies
   type(physics_buffer_desc), pointer       :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)             ! Land fraction
   integer, intent(in) :: ideep(pcols)
   
   ! Local variables

   integer :: i,k
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: ixnumice, ixnumliq      ! constituent indices for cloud liquid and ice number concentration. 
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns

   real(r8) :: ftem(pcols,pver)       ! Temporary workspace for outfld variables
   real(r8) :: num_tem                ! Temporary holder for number tendency
   
   real(r8) :: q_delta_adv(pcols,pver), q_delta_auto(pcols,pver), q_delta_sed(pcols,pver), t_delta_adv(pcols,pver), t_delta_auto(pcols,pver), t_delta_sed(pcols,pver), t_delta_rad(pcols,pver), t_flux_adv(pcols,pver)
   real(r8) :: precflx(pcols,pver), snowflx(pcols,pver)

   integer :: lengath
   real(r8) :: ptends(pcols,pver),ptendqv(pcols,pver),ptendqi(pcols,pver),ptendql(pcols,pver)
   real(r8) :: q_delta_adv_(pcols,pver), q_delta_auto_(pcols,pver), q_delta_sed_(pcols,pver), t_delta_adv_(pcols,pver), t_delta_auto_(pcols,pver), t_delta_sed_(pcols,pver), t_delta_rad_(pcols,pver), t_flux_adv_(pcols,pver), precflx_(pcols,pver), snowflx_(pcols,pver)
   real(r8) :: precsfc(pcols),snowsfc(pcols)
   
   real(r8), pointer, dimension(:,:) :: dp_cldice, dp_cldliq, dp_flxprc, dp_flxsnw
   logical  :: lq(pcnst)

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! total snow
   
   real(r8) :: yog_precsfc(pcols)   ! precipitation at surface
   real(r8) :: yog_snowsfc(pcols)   ! snow at surface
   
   lchnk = state%lchnk
   ncol  = state%ncol

   !write(iulog,*)'ideep',ideep(:ncol)
   
   ! fetch prec_dp from the physics buffer
   ! Note that prec_dp is initialised in the physics buffer and zeroed for deep_scheme='off'
   prec_dp_idx     = pbuf_get_index('PREC_DP')
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   snow_dp_idx     = pbuf_get_index('SNOW_DP')
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )   

   dp_flxprc_idx = pbuf_get_index('DP_FLXPRC')
   call pbuf_get_field(pbuf, dp_flxprc_idx, dp_flxprc )
   dp_flxsnw_idx = pbuf_get_index('DP_FLXSNW')
   call pbuf_get_field(pbuf, dp_flxsnw_idx, dp_flxsnw)
   dp_cldliq_idx = pbuf_get_index('DP_CLDLIQ')
   call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq )
   dp_cldice_idx = pbuf_get_index('DP_CLDICE')
   call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice )

   dp_cldliq(:ncol,:) = 0._r8
   dp_cldice(:ncol,:) = 0._r8
  
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
      
   lq(:) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'yogNN', ls=.true., lq=lq(:))! initialize ptend type for YOG

   call nn_convection_flux_CAM(state%pmid(:,pver:1:-1), state%pint(:,pverp:1:-1), state%ps, &
                               state%t(:,pver:1:-1), state%q(:,pver:1:-1,1), &
                               state%q(:,pver:1:-1,ixcldliq), state%q(:,pver:1:-1,ixcldice), &
                               state%phis, state%zm(:,pver:1:-1), &
                               landfrac, state%ps, &
                               cpair, &
                               ztodt, &
                               ncol, &
                               yog_precsfc, yog_snowsfc, &
                               ptend%q(:,pver:1:-1,ixcldice), ptend%q(:,pver:1:-1,1), &
                               ptend%q(:,pver:1:-1,ixcldliq), ptend%s(:,pver:1:-1), &
                               q_delta_adv(:,pver:1:-1), q_delta_auto(:,pver:1:-1), &
                               q_delta_sed(:,pver:1:-1), t_delta_adv(:,pver:1:-1), &
                               t_delta_auto(:,pver:1:-1), t_delta_sed(:,pver:1:-1), &
                               t_delta_rad(:,pver:1:-1), t_flux_adv(:,pver:1:-1),&
                               precflx(:,pver:1:-1), snowflx(:,pver:1:-1))
   
   ! XL XL XL
   lengath = count(ideep > 0)
   ptends(:,:) = 0._r8
   ptendqv(:,:) = 0._r8
   ptendqi(:,:) = 0._r8
   ptendql(:,:) = 0._r8
   precsfc(:) = 0._r8
   snowsfc(:) = 0._r8

   precflx_(:,:) = 0._r8
   snowflx_(:,:) = 0._r8
   
   q_delta_adv_(:,:) = 0._r8
   q_delta_auto_(:,:) = 0._r8
   q_delta_sed_(:,:) = 0._r8
   t_delta_adv_(:,:) = 0._r8
   t_delta_auto_(:,:) = 0._r8
   t_delta_sed_(:,:) = 0._r8
   t_delta_rad_(:,:) = 0._r8
   t_flux_adv_(:,:) = 0._r8
      
   do i = 1,lengath
      do k=1,pver
      ptends(ideep(i),k) = ptend%s(ideep(i),k)
      ptendqv(ideep(i),k) = ptend%q(ideep(i),k,1)
      ptendqi(ideep(i),k) = ptend%q(ideep(i),k,ixcldice)
      ptendql(ideep(i),k) = ptend%q(ideep(i),k,ixcldliq)      
      precsfc(ideep(i)) = yog_precsfc(ideep(i))
      snowsfc(ideep(i)) = yog_snowsfc(ideep(i))

      q_delta_adv_(ideep(i),k) = q_delta_adv(ideep(i),k)
      q_delta_auto_(ideep(i),k) = q_delta_auto(ideep(i),k)
      q_delta_sed_(ideep(i),k) = q_delta_sed(ideep(i),k)
      t_delta_adv_(ideep(i),k) = t_delta_adv(ideep(i),k)
      t_delta_auto_(ideep(i),k) = t_delta_auto(ideep(i),k)
      t_delta_sed_(ideep(i),k) = t_delta_sed(ideep(i),k)
      t_delta_rad_(ideep(i),k) = t_delta_rad(ideep(i),k)
      t_flux_adv_(ideep(i),k) = t_flux_adv(ideep(i),k)

      precflx_(ideep(i),k) = precflx(ideep(i),k)
      snowflx_(ideep(i),k) = snowflx(ideep(i),k)
      
      end do
   end do
   
   ftem(:ncol,:pver) = ptends(:ncol,:pver)/cpair
   call outfld('YOGDT   ',ftem ,pcols   ,lchnk   )
   call outfld('YOGDQ   ',ptendqv ,pcols  ,lchnk   )
   call outfld('YOGDICE ',ptendqi ,pcols   ,lchnk   )
   call outfld('YOGDLIQ ',ptendql ,pcols   ,lchnk   )
   call outfld('PREC_YOG',precsfc ,pcols   ,lchnk   )
   call outfld('SNOW_YOG',snowsfc ,pcols   ,lchnk   )   

   call outfld('YOG_DQ_ADV   ',q_delta_adv_ ,pcols   ,lchnk   )
   call outfld('YOG_DQ_AUTO  ',q_delta_auto_,pcols   ,lchnk   )
   call outfld('YOG_DQ_SED   ',q_delta_sed_ ,pcols   ,lchnk   )
   call outfld('YOG_DT_ADV   ',t_delta_adv_ ,pcols   ,lchnk   )
   call outfld('YOG_DT_AUTO  ',t_delta_auto_,pcols   ,lchnk   )        
   call outfld('YOG_DT_SED   ',t_delta_sed_ ,pcols   ,lchnk   )
   call outfld('YOG_DT_RAD   ',t_delta_rad_,pcols   ,lchnk   )

   call outfld('YOG_DT_ADVFLX   ', t_flux_adv_,pcols   ,lchnk   )           


   do k = 1, pver
      do i = 1, ncol

      ptend%s(i,k) = ptends(i,k)
      ptend%q(i,k,1) = ptendqv(i,k)
      ptend%q(i,k,ixcldice) = ptendqi(i,k)
      ptend%q(i,k,ixcldliq) = ptendql(i,k)
      
      end do
   end do
   
   ! Add prec_yog to prec_dp for passing to physpkg and coupler
   ! Note that prec_dp is initialised in the physics buffer and zeroed for deep_scheme='off'
   prec(:ncol) = prec(:ncol) + precsfc(:ncol)
   snow(:ncol) = snow(:ncol) + snowsfc(:ncol)

   dp_flxprc(:ncol,:pver) = precflx_(:ncol,:pver)
   dp_flxsnw(:ncol,:) = snowflx_(:ncol,:pver)
   
   ! Update the number concentration tendencies for liquid and ice species for all cells
   ! Match calculations in the `clubb_tend_cam()` subroutine
   ! YOG could produce negative number tendency, so we ensure it doesn't reduce total
   ! number concentration below 0.0 - note that a check on this is also made in physics_update()
   ! XL XL XL
   call cnst_get_ind('NUMLIQ', ixnumliq)
   call cnst_get_ind('NUMICE', ixnumice)

   do k = 1, pver
      do i = 1, ncol
         ! Liquid
         if (ptend%q(i,k,ixcldliq) .lt. 0) then
            if (state%q(i,k,ixcldliq) .gt. 0) then
               ! change number in proportion to current number
               num_tem = ptend%q(i,k,ixcldliq) / state%q(i,k,ixcldliq) * state%q(i,k,ixnumliq)
            else
               num_tem = 0.0
            endif
            ptend%q(i,k,ixnumliq) = - min(-num_tem, state%q(i,k,ixnumliq)/ztodt)
         else
            num_tem = 3. * ptend%q(i,k,ixcldliq) / (4.0*3.14* 8.0e-6**3*997.0)
            ptend%q(i,k,ixnumliq) = num_tem
         endif
         ! Ice
         if (ptend%q(i,k,ixcldice) .lt. 0) then
            if (state%q(i,k,ixcldice) .gt. 0) then
               ! change number in proportion to current number
               num_tem = ptend%q(i,k,ixcldice) / state%q(i,k,ixcldice) * state%q(i,k,ixnumice)
            else
               num_tem = 0.0
            endif
            ptend%q(i,k,ixnumice) = - min(-num_tem, state%q(i,k,ixnumice)/ztodt)
         else
            num_tem = 3. * ptend%q(i,k,ixcldice) / (4.0*3.14*25.0e-6**3*500.0) * 0.5  ! XL XL XL
            ptend%q(i,k,ixnumice) = num_tem
         endif
      end do
   end do

   call outfld('YOGDNUMLIQ ',ptend%q(1,1,ixnumliq) ,pcols   ,lchnk   )
   call outfld('YOGDNUMICE ',ptend%q(1,1,ixnumice) ,pcols   ,lchnk   )
   ! XL XL XL
   
end subroutine yog_tend

!=========================================================================================

end module yog_intr
