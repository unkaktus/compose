!***********************************************************************
!***********************************************************************
! modules
!***********************************************************************
!***********************************************************************
MODULE eos_tables
! Stefan Typel for the CompOSE core team, version 1.11, 2017/05/22

! parameters
! dim_reg:   number of standard thermodynamic quantities in eos.thermo (7)
! dim_qtyt:  regular number of quantities in array eos_thermo (23)
! dim_qtye:  maximum number of error quantities (8)

implicit none

integer, parameter :: dim_reg=7,dim_qtyt=23,dim_qtye=8,dim_df=10

integer :: nall_max,nadd_max,np_max,nq_max,nm_max

integer, allocatable :: idxp_compo(:,:,:,:),idxq_compo(:,:,:,:),idx_mic(:,:,:,:),&
     idx_compo_p(:),idx_compo_q(:),idx_micro(:)

double precision :: arg(3),eos_thermo(dim_qtyt),eos_thermo_err(dim_qtye),&
     eos_df(dim_df)

double precision, allocatable :: tab_para(:,:),&
     tab_thermo(:,:,:,:),eos_thermo_add(:),&
     tabp_compo(:,:,:,:),tabq_compo(:,:,:,:),tab_mic(:,:,:,:),&
     eos_compo_p(:),eos_compo_q(:,:),eos_micro(:),eos_q(:)

! neutron mass in MeV/c^{2}, standard value
double precision :: m_n =  939.565379d00

! proton mass in MeV/c^{2}, standard value
double precision :: m_p =  938.272046d00

END MODULE eos_tables
!***********************************************************************
MODULE compose_internal
! Stefan Typel for the CompOSE core team, version 1.14, 2017/05/22
USE eos_tables
implicit none

integer, parameter :: dim_err=100,dim_ip=25,dim_iq=10,&
     dim_im=10,dim_ia=5,nerr_max=8

integer :: eos_dim,idx_ex(0:3),idx_ipl(3),&
     dim_a,dim_idx(3),min_idx(4),&
     idx_arg2(-4:5,-4:5,0:2),idx_arg1(-4:5,0:2),imap(3),jmap(4),&
     ipl_rule(0:1,3,-8:8),&
     error_msg(0:dim_err),&
     n_qty,n_add,idx_qty(dim_qtyt),idx_a(dim_ia),&
     n_df,idx_df(dim_df),&
     n_p,n_q,n_m,n_err,idx_err(dim_qtye),iout,&
     incl_l,inbyq,dim_p(4),irpl

integer, allocatable :: idx_arg(:,:,:,:),idx_argx(:,:,:,:),idx_add(:),&
     idx_thermo(:),idx_p(:),idx_q(:),idx_m(:)

double precision :: d1x(-4:5,-4:5,-4:4),d2x(-4:5,-4:5,-4:4),&
     d1y(-4:5,-4:5,-4:4),d2y(-4:5,-4:5,-4:4),dx,dy,dx2,dy2,&
     df(0:3,0:3,-4:5,-4:5),fc(0:5,0:5),val_rpl,&
     para_min(4),para_max(4),arg2(4)
 
double precision, allocatable :: r1d(:,:,:,:),r2d(:,:,:,:),v_thermo(:,:)

END MODULE compose_internal
!***********************************************************************
