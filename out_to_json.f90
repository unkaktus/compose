!***********************************************************************
!***********************************************************************
! Module to export the column attribution in eos.table in a json
! format for web use in the CompOSE project.
! CompOSE core team, version 1.0, 2017/12/07
!***********************************************************************
!***********************************************************************


module  m_out_to_json

 implicit none

 private

 public :: out_to_json_write

contains

 subroutine out_to_json_write(iwr,irpl,&
   &                          idx_qty, idx_add, idx_df, idx_p, idx_q, idx_m, idx_err,&
   &                            n_qty,   n_add,    n_df,  n_p,   n_q,   n_m,   n_err)
  integer, intent(in) :: iwr, n_qty, n_add, n_df, n_p, n_q, n_m, n_err, irpl
  integer, intent(in),dimension(:) :: idx_qty, idx_add, idx_df, idx_p, idx_q, idx_m, idx_err
  integer :: i1,i2,ifile = 10000
  character(11) :: f1 = '(a,i0,a)'
  character(16) :: f2 = '(a,i0,a,i6,a)'

  open(unit=ifile,file='eos.info.json',status='replace')

  if (iwr /= 1) return
  i2 = 0
  write(ifile,*) '{ '
  write(ifile,*) '  "info" : "The columns of the file eos.table contains the following quantities",'
  write(ifile,*) '  "columns" : {'
  i2 = i2+1
  write(ifile,f1) '"',i2,'":{"title": "temperature T", "unit": "MeV" ,"symbol":"T"},'
  i2 = i2+1
  write(ifile,f1) '"',i2,'":{"title": "baryon number density n_b", "unit": "fm^-3","symbol":"nb"},'
  i2 = i2+1
  write(ifile,f1) '"',i2,'":{"title": "hadronic charge fraction Y_q", "unit":"","symbol":"Yq"},'
  if (irpl > 0) then
    i2 = i2+1
    write(ifile,f1) '"',i2,'":{"title": "magnetic field strength B", "unit":"G","symbol":"B"},'
  endif
  do i1=1,n_qty,1
    i2 = i2+1
    select case(idx_qty(i1))
    case(1)
      write(ifile,f1) '"',i2,'":{"title": "pressure p", "unit":"MeV fm^-3"},'
    case(2)
      write(ifile,f1) '"',i2,'":{"title": "entropy per baryon S", "unit":""},'
    case(3)
      write(ifile,f1) '"',i2,'":{"title": "shifted baryon chemical potential mu_b-m_n", "unit":"MeV"},'
    case(4)
      write(ifile,f1) '"',i2,'":{"title": "charge chemical potential mu_q", "unit":"MeV"},'
    case(5)
      write(ifile,f1) '"',i2,'":{"title": "lepton chemical potential mu_l", "unit":"MeV"},'
    case(6)
      write(ifile,f1) '"',i2,'":{"title": "scaled free energy per baryon F/m_n-1", "unit":""},'
    case(7)
      write(ifile,f1) '"',i2,'":{"title": "scaled internal energy per baryon E/m_n-1", "unit":""},'
    case(8)
      write(ifile,f1) '"',i2,'":{"title": "scaled enthalpy energy per baryon H/m_n-1", "unit":""},'
    case(9)
      write(ifile,f1) '"',i2,'":{"title": "scaled free enthalpy per baryon G/m_n-1", "unit":""},'
    case(10)
      write(ifile,f1) '"',i2,'":{"title": "derivative dp/dn_b|E", "unit":"MeV"},'
    case(11)
      write(ifile,f1) '"',i2,'":{"title": "derivative p/dE|n_b", "unit":"fm^-3"},'
    case(12)
      write(ifile,f1) '"',i2,'":{"title": "square of speed of sound (c_s)^2 ", "unit":""},'
    case(13)
      write(ifile,f1) '"',i2,'":{"title": "specific heat capacity at constant volume c_V", "unit":""},'
    case(14)
      write(ifile,f1) '"',i2,'":{"title": "specific heat capacity at constant pressure c_p", "unit":""},'
    case(15)
      write(ifile,f1) '"',i2,'":{"title": "adiabatic index Gamma", "unit":""},'
    case(16)
      write(ifile,f1) '"',i2,'":{"title": "expansion coefficient at constant pressure alpha_p ", "unit":"MeV^-1"},'
    case(17)
      write(ifile,f1) '"',i2,'":{"title": "tension coefficient at constant volume beta_V", "unit":"fm^-3"},'
    case(18)
      write(ifile,f1) '"',i2,'":{"title": "isothermal compressibility kappa_T", "unit":"MeV^-1 fm^3"},'
    case(19)
      write(ifile,f1) '"',i2,'":{"title": "isentropic compressibility kappa_S", "unit":"MeV^-1 fm^3"},'
    case(20)
      write(ifile,f1) '"',i2,'":{"title": "free energy per baryon F", "unit":"MeV"},'
    case(21)
      write(ifile,f1) '"',i2,'":{"title": "internal energy per baryon E", "unit":"MeV"},'
    case(22)
      write(ifile,f1) '"',i2,'":{"title": "enthalpy per baryon H", "unit":"MeV"},'
    case(23)
      write(ifile,f1) '"',i2,'":{"title": "free enthalpy per baryon G", "unit":"MeV"},'
    end select
  end do
  if (n_add > 0) then
    do i1=1,n_add,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "additional thermodynamic quantity", "index":',idx_add(i1),'},'
    end do
  end if
  !2017/05/22
  if (n_df > 0) then
    do i1=1,n_df,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "derivative quantity", "index":',idx_df(i1),'},'
    end do
  end if
  if (n_p > 0) then
    do i1=1,n_p,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "number fraction Y of particle with index", "index":',idx_p(i1),'},'
    end do
  end if
  if (n_q > 0) then
    do i1=1,n_q,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "total number fraction Y of particle set with index", "index":',idx_q(i1),'},'
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "average mass number A_av of particle set with index", "index":',idx_q(i1),'},'
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "average proton number Z_av of particle set with index", "index":',idx_q(i1),'},'
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "average neutron number N_av of particle set with index", "index":',idx_q(i1),'},'
    end do
  end if
  if (n_m > 0) then
    do i1=1,n_m,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "microscopic quantity with index", "index":',idx_m(i1),'},'
    end do
  end if
  if (n_err > 0) then
    do i1=1,n_err,1
      i2 = i2+1
      write(ifile,f2) '"',i2,'":{"title": "error quantity with index", "index":',idx_err(i1),'},'
    end do
  end if
  write(ifile,*) '   "-1": {"title": "closing element", "index": -1}'
  write(ifile,*)'  }'
  write(ifile,*)'}'

  close(unit=ifile)

 end subroutine out_to_json_write



end module m_out_to_json
