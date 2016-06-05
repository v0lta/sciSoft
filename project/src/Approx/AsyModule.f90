!Asymptotic expansion computation module. This module contains translated
!matlab code from http://nines.cs.kuleuven.be/software/JACOBI/ and toledo.
module AsyModule
  use UtilModule
  implicit none

contains

  subroutine asyLens(n,nrT,z,alpha,beta,Dinf,Uright,Uleft,ipi)
    integer, intent(in) :: n,nrT
    complex(kind = PRECISION), intent(in) :: z
    real(kind = PRECISION), intent(in) :: alpha,beta,Dinf
    complex(kind = PRECISION), intent(in), dimension(:,:,:,:) :: Uright,Uleft
    complex(kind = PRECISION) ::ipi, nmz = 0
    complex(kind = PRECISION), dimension(2,2) :: RI = 0
    integer :: i,k,m
    complex(kind = PRECISION), dimension(2) :: tmp1
    complex(kind = PRECISION), dimension(2) :: tmp2
    logical, parameter :: debug = .false.

    if (debug .eqv. .true.) then
        print *, 'n: ', n, 'nrT: ', nrt, 'z: ', z, 'alpha: ', alpha, 'beta: ', &
        beta, 'Dinf:', Dinf
        print *,(/ ('-',i=1,60) /)
    end if


    if ( (nrT .eq. 1) .or. (size(Uright) .eq. 0)) then
      nmz = 1
    else
      nmz = sqrt(real(1.0+2.0*cmplx(0.0,1.0)*Dinf**2*sum(reshape(Uright(2,1,1:(nrT-1),1) + &
            Uleft(2,1,1:(nrT-1),1), (/ (nrT-1), 1 /)) &
            /(n+1)**reshape((/ (i,i=1,(nrT-1)) /), (/ (nrT-1), 1 /) ))));
    end if
    if (debug .eqv. .true.) then
      print *, 'nmz:'
      print *, nmz
      print *, 'Uright:', Uright(2,1,1:(nrT-1),1)
      print *, 'Uleft:',   Uleft(2,1,1:(nrT-1),1)
    end if

    call eye(cmplxMatrix=RI)
    do k = 1,(nrT-1)
      do m = 1,nint(k/2.0)
        RI = RI + (Uright(:,:,k,m)/(z-1)**m + Uleft(:,:,k,m)/(z+1)**m)/n**k;
      end do
    end do

    if (debug .eqv. .true.) then
      print *, 'RI:'
      call printRows(cmplxMatrix=RI)
    end if


    tmp1 = matmul((/ 2**(1.0/2.0-n) /( (1.0+z)**(1.0/4.0+alpha/2.0)*(1.0-z) &
           **(1.0/4.0+beta/2.0) ), cmplx(0.0,0.0) /),RI)
                                       !ifort does not allow complx acos input.     
    tmp2 = (/ cmplx(Dinf,0.0)  *  cos((n+1.0/2.0)*acos(real(z)) + (alpha+beta)/2.0 - &
           (alpha * pi)/2.0  -pi/4.0), &
           cmplx(0.0,-1/Dinf)*cos((n-1.0/2.0)*acos(real(z)) + (alpha+beta)/2.0 - &
           (alpha * pi)/2.0-pi/4.0) /)
    ipi = dot_product(conjg(tmp1),tmp2)*nmz

    if (debug .eqv. .true.) then
        print *, 'tmp1:', tmp1
        print *, 'tmp2:', tmp2
        print *, dot_product(CONJG(tmp1),tmp2)
        
        print *, 'ipi: ', ipi
    end if
  end subroutine asyLens !----------------------------------------------------------------

  subroutine compUExplicit(alpha, beta, Dinf, c, d, Uright, Uleft)
    real(kind = PRECISION), intent(in) :: alpha
    real(kind = PRECISION), intent(in) :: beta
    real(kind = PRECISION), intent(in) :: Dinf
    complex(kind = PRECISION), dimension(:), intent(in) :: c
    complex(kind = PRECISION), dimension(:), intent(in) :: d
    complex(kind = PRECISION), intent(out), dimension(:,:,:,:) :: Uleft
    complex(kind = PRECISION), intent(out), dimension(:,:,:,:) :: Uright
    complex(kind = PRECISION) :: preFactor
    complex(kind = PRECISION), dimension(2,2) :: tmpStore
    logical, parameter :: debug = .false.

    Uright = 0
    Uleft = 0
    
    if (debug .eqv. .true.) then
        print *, '----------------------- Uexp Input'
        print *, 'alpha: ', alpha, 'beta: ', beta , 'Dinf: ', Dinf, 'c: ', c , 'd: ', d
        print *, '----------------------------------'
    end if

    preFactor = (4.0*alpha**2.0-1.0)*0.0625
    tmpStore =matmul(matmul(transpose(reshape((/ cmplx(Dinf,0.0) ,cmplx(0.0,0.0),  &
              cmplx(0.0,0.0), cmplx(Dinf,0.0)**(-1) /), (/ 2,2 /))),               &
              transpose(reshape((/ cmplx(-1.0,0.0) ,cmplx(0.0,1.0),cmplx(0.0,1.0), & 
              cmplx(1.0,0.0) /), (/ 2,2 /)))),                                     &
              transpose(reshape((/ cmplx(Dinf,0.0)**(-1),cmplx(0.0,0.0),           &
              cmplx(0.0,0.0),cmplx(Dinf,0.0) /), (/2,2 /))))
    Uright(:,:,1,1) = preFactor * tmpStore
    
    if (debug .eqv. .true.) then
        print *, 'Uright(:,:,1,1)::prefactor: ', preFactor
        !call printRows(cmplxMatrix=tmpStore)
        print *, 'Uright(:,:,1,1)'
        call printRows(cmplxMatrix=Uright(:,:,1,1))
    end if

    preFactor = (4.0*beta**2.0-1.0)*0.0625
    tmpStore = matmul(matmul(transpose(reshape((/ cmplx(Dinf,0.0),cmplx(0.0,0.0), &
               cmplx(0.0,0.0),cmplx(Dinf,0.0)**(-1.0) /), (/ 2,2 /))) ,           &
               transpose(reshape((/ cmplx(1.0,0.0),cmplx(0.0,1.0),cmplx(0.0,1.0), &
               cmplx(-1.0,0.0) /),(/ 2,2 /)))),                                   &
               transpose(reshape((/ cmplx(Dinf,0.0)**(-1.0),cmplx(0.0,0.0),       &
               cmplx(0.0,0.0),cmplx(Dinf,0.0) /), (/2,2 /))))
    Uleft(:,:,1,1) = preFactor * tmpStore
    
    if (debug .eqv. .true.) then
        print *, 'Uleft(:,:,1,1)::prefactor: ', preFactor
        !call printRows(cmplxMatrix=tmpStore)
        print *, 'Uleft(:,:,1,1)'
        call printRows(cmplxMatrix=Uleft(:,:,1,1))
    end if


    preFactor = (4.0*alpha**2.0-1.0)*0.00390625
    tmpStore =  matmul(matmul(transpose(reshape((/ cmplx(Dinf,0.0),cmplx(0.0,0.0),    & 
                cmplx(0.0,0.0),cmplx(Dinf,0.0)**(-1.0) /), (/ 2,2 /))),               &
                transpose(reshape((/ 8.0*alpha+8.0*beta+8.0*c(0+1)-4*beta**2.0+1.0,   &
                cmplx(0.0,1.0)*(-8.0*alpha-8.0*beta-8.0*c(0+1)+4.0*alpha**2.0 +4.0    &
                *beta**2.0-10.0),                                                     &
                cmplx(0.0,1.0)*(-8.0*alpha-8.0*beta-8.0*c(0+1)-4.0*alpha**2.0 -4.0    &
                *beta**2.0+10.0),                                                     &
                -8.0*alpha-8.0*beta-8.0*c(0+1)-4.0*beta**2.0+1.0 /), (/ 2,2 /)))),    &
                transpose(reshape((/ cmplx(Dinf,0.0)**(-1),cmplx(0.0,0.0),            &
                cmplx(0.0,0.0),cmplx(Dinf,0.0) /), (/ 2,2 /))));
    Uright(:,:,2,1) = preFactor * tmpStore                                                            
    
    if (debug .eqv. .true.) then
        print *, 'Uright(:,:,2,1)::prefactor: ', preFactor
        !call printRows(cmplxMatrix=tmpStore)
        print *, 'Uright(:,:,2,1):'
        call printRows(cmplxMatrix=Uright(:,:,2,1))
    end if


    preFactor = (4.0*beta**2.0-1.0)*0.00390625
    tmpStore = matmul(matmul( transpose(reshape((/ cmplx(Dinf,0.0),cmplx(0.0,0.0),    & 
               cmplx(0.0,0.0),cmplx(Dinf,0.0)**(-1.0) /), (/ 2,2 /))),                &
               transpose(reshape((/ cmplx(-1.0,0.0)*(+8.0*beta + 8.0*alpha-8.0*d(0+1) &
               - 4*alpha**2.0  + 1.0),                                                &
               cmplx(0.0,1.0) *(-8.0*beta - 8.0*alpha+8.0*d(0+1) + 4.0*beta**2.0 +    & 
               4.0*alpha**2.0-10.0),                                                  &
               cmplx(0.0,1.0) *(-8.0*beta - 8.0*alpha+8.0*d(0+1) - 4.0*beta**2.0 -    &
               4.0*alpha**2.0+10.0),                                                  &
               cmplx(-1.0,0.0)*(-8.0*beta - 8.0*alpha+8.0*d(0+1) - 4.0*alpha**2.0+    &
               1.0) /), (/ 2,2 /)))),                                                 &
               transpose(reshape((/ cmplx(Dinf,0.0)**(-1.0),cmplx(0.0,0.0),           &
               cmplx(0.0,0.0), cmplx(Dinf,0.0) /), (/ 2,2 /))))
    Uleft(:,:,2,1) = preFactor * tmpStore  
    
    if (debug .eqv. .true.) then
       print *, 'Uleft(:,:,2,1)::prefactor: ', preFactor
       !call printRows(cmplxMatrix=tmpStore)
       print *, 'Uleft(:,:,2,1):'
       call printRows(cmplxMatrix=Uleft(:,:,2,1))
    end if                                                        
    

    Uright(:,:,3,1) = transpose(reshape((/ 0.0078125*d(0+1)*alpha**2-0.0078125*c(0+1)*   &
      beta**2+0.0625*alpha**3*beta**2+0.0625*alpha**2*beta**3-0.0078125*alpha**4*beta**  &
      2-0.0078125*alpha**2*beta**4-0.015625*alpha*beta**2+0.0078125*d(0+1)*beta**2-0.125 &
      *alpha**3*beta-0.015625*alpha**2*beta-0.0625*alpha**2*c(0+1)**2-0.125*c(0+1)       &
      *alpha**3+0.03125*c(0+1)*alpha+0.03125*beta*c(0+1)-0.0546875*alpha**2*beta**2      &
      -0.0078125*c(0+1)*alpha**2+0.03125*alpha*beta+0.03125*c(0+1)*alpha**2*beta**2-     &
      0.03125*d(0+1)*alpha**2*beta**2-0.125*c(0+1)*alpha**2*beta+0.001953125*c(0+1)-     &
      0.001953125*d(0+1)+0.000244140625+0.001953125*beta**4-0.015625*beta**3-            &
      0.015625*alpha**3+0.015625*c(0+1)**2-0.060546875*alpha**4+0.00390625*alpha         &
      +0.00390625*beta+0.01416015625*alpha**2+0.01416015625*beta**2,                     &
      cmplx(0.0,(1.0/6144.0))*Dinf**2*(-132-48*d(0+1)*alpha**2+48*c(0+1)*beta**2         &
      -384*alpha**3*beta**2-384*alpha**2*beta**3+144*alpha**4*beta**                     &
      2+48*alpha**2*beta**4+96*alpha*beta**2-48*d(0+1)*beta**2+768*alpha**3*beta+        &
      1056*alpha**2*beta+384*alpha**2*c(0+1)**2+768*c(0+1)*alpha**3-192*c(0+1)*alpha-    &
      192*beta*c(0+1)+96*alpha**2*beta**2+1008*c(0+1)*alpha**2-384*alpha**4*beta-384*    &
      c(0+1)*alpha**4-192*alpha*beta-192*c(0+1)*alpha**2*beta**2+192*d(0+1)*alpha**2*    &
      beta**2+768*c(0+1)*alpha**2*beta-228*c(0+1)+12*d(0+1)-384*alpha**5+64*alpha**6     &
      -12*beta**4+96*beta**3+1056*alpha**3-96*c(0+1)**2-20*alpha**4-240*alpha-240*beta   &
      +529*alpha**2-33*beta**2),                                                         &                                          
      cmplx(0.0,(-1.0/6144.0))*(132 -48 *d(0+1)*alpha**2+48*c(0+1) *beta**2-384*         &
      alpha**3*beta**2-384*alpha**2*beta**3-144*alpha**4*beta**2-48*alpha**2*beta**4+    &
      96*alpha*beta**2-48*d(0+1)*beta**2-768*alpha**3*beta+1056*alpha**2*beta-384*       &
      alpha**2*c(0+1)**2-768*c(0+1)*alpha**3+192*c(0+1)*alpha+192*beta*c(0+1)-96*        &
      alpha**2*beta**2+1008*c(0+1)* alpha**2-384*alpha**4*beta-384*c(0+1)*alpha**        &
      4+192*alpha*beta-192*c(0+1)*alpha**2*beta**2+192*d(0+1)*alpha**2*beta**2-768*      &
      c(0+1)*alpha**2*beta-228*c(0+1)+12*d(0+1)-384*alpha**5-64*alpha**6+12*beta**4      &
      +96*beta**3+1056*alpha**3+96*c(0+1)**2+20*alpha**4-240*alpha-240*beta-529*alpha**2 &
      +33*beta**2)/Dinf**2,                                                              &
      0.0078125*d(0+1)*alpha**2-0.0078125*c(0+1)*beta**2+0.0625*alpha**3*beta**2         &
      +0.0625*alpha**2*beta**3+0.0078125*alpha**4*beta**2+0.0078125*alpha**2*beta**4     &
      -0.015625*alpha*beta**2+0.0078125*d(0+1)*beta**2+0.125*alpha**3*beta-0.015625*     &
      alpha**2*beta+0.0625*alpha**2*c(0+1)**2+0.125*c(0+1)*alpha**3-0.03125*c(0+1)       &
      *alpha-0.03125*beta*c(0+1)+0.0546875*alpha**2*beta**2-0.0078125*c(0+1)*alpha**2    &
      -0.03125*alpha*beta+0.03125*c(0+1)*alpha**2*beta**2-0.03125*d(0+1)*alpha**2        &
      *beta**2+0.125*c(0+1)*alpha**2*beta+0.001953125*c(0+1)-0.001953125*                &
      d(0+1)-0.000244140625-0.001953125*beta**4-0.015625*beta**3-0.015625*               &
      alpha**3-0.015625*c(0+1)**2+0.060546875*alpha**4+0.00390625*alpha+0.00390625*      &
      beta-0.01416015625*alpha**2-0.01416015625*beta**2 /) , (/ 2,2 /))) 

    if (debug .eqv. .true.) then
       print *, 'Uright(:,:,3,1):'
       call printRows(cmplxMatrix=Uright(:,:,3,1))
    end if        


    Uright(:,:,3,2) = transpose(reshape((/ cmplx((-1.0/192.0),0.0)*alpha**6             &
      +(35.0/768.0)*alpha**4-(259.0/3072.0)*alpha**2+(75.0/4096.0),(1.0/192.0)*         &
      alpha**6-(35.0/768.0)*alpha**4+(259.0/3072.0)*alpha**2-(75.0/4096.0)              &
      *cmplx(0.0,1.0)*Dinf**2.0, cmplx(0.0,(1.0/12288.0))*(64.0*alpha**6.0-560.0*       &
      alpha**4.0+1036.0*alpha**2.0-225.0)/Dinf**2.0, cmplx((1.0/192.0),0.0)*            &
      alpha**6.0-(35.0/768.0)*alpha**4.0+(259.0/3072.0)*                                &
      alpha**2.0-(75.0/4096.0) /), (/ 2,2 /))) !ok
      
    if (debug .eqv. .true.) then
       print *, 'Uright(:,:,3,2):'
       call printRows(cmplxMatrix=Uright(:,:,3,2))
    end if  

    Uleft(:,:,3,1) = transpose(reshape((/ cmplx(0.03125,0.0)*beta*d(0+1)                &
      +0.03125*alpha*d(0+1)+0.0625*beta**2*d(0+1)**2-(1.0/8.0)                          &
      *d(0+1)*beta**3+0.125*alpha*beta**3-0.0078125*d(0+1)*alpha**2+0.0078125           &
      *c(0+1)*beta**2-0.0625*alpha**3*beta**2-0.0625*alpha**2*beta**3+0.0078125*        &
      alpha**4*beta**2+0.0078125*alpha**2*beta**4+(1.0/64.0)*alpha*beta**2.0            &
      -0.0078125*d(0+1)*beta**2+0.015625*alpha**2*beta+7/128*alpha**2*beta**2           &
      +(1.0/128.0)*c(0+1)*alpha**2.0-0.03125*alpha*beta-0.03125*c(0+1)*alpha**2.0*      &
      beta**2+0.03125*d(0+1)*alpha**2.0*beta**2.0-0.125*d(0+1)*alpha*                   &
      beta**2.0-0.001953125*c(0+1)+0.001953125*d(0+1)-0.000244140625-                   &
      0.015625*d(0+1)**2.0+0.060546875* beta**4.0+0.015625*beta**3.0+0.015625*          &
      alpha**3.0-0.001953125*alpha**4-0.00390625*alpha-0.00390625*beta                  &
      -0.01416015625* alpha**2.0-0.01416015625*beta**2.0,                               &
      cmplx(0.0,1.0/6144.0)*Dinf**2*(-132-48*d(0+1)*alpha**2+48*c(0+1)*beta**2-384*     &
      alpha**3*beta**2-384*alpha**2*beta**3+144*alpha**4*beta**2+48*alpha**2*           &
      beta**4+96*alpha*beta**2-48*d(0+1)*beta**2+768*alpha**3*beta+1056*alpha**2*       &
      beta+384*alpha**2*c(0+1)**2+768*c(0+1)*alpha**3-192*c(0+1)*alpha-192*beta*        &
      c(0+1)+96*alpha**2*beta**2+1008*c(0+1)*alpha**2-384*alpha**4*beta-384*c(0+1)*     &
      alpha**4-192*alpha*beta-192*c(0+1)*alpha**2*beta**2+192*d(0+1)*alpha**2*beta**2   &
      +768*c(0+1)*alpha**2*beta-228*c(0+1)+12*d(0+1)-384*alpha**5+64*alpha**6-12*       &
      beta**4+96*beta**3+1056*alpha**3-96*c(0+1)**2-20*alpha**4-240*alpha-240*          &
      beta+529*alpha**2-33*beta**2) ,                                                   &
      cmplx(0.0,(1.0/6144.0))*(-132.0+192.0*beta*d(0+1)+192.0*alpha*d(0+1)+384.0*       &
      beta**2.0*d(0+1)**2.0-768.0*d(0+1)*beta**3.0+768.0*alpha*beta**3.0-384.0*d(0+1)*  &
      beta**4.0+384.0*alpha*beta**4.0+48.0*d(0+1)*alpha**2.0-48.0*c(0+1)*beta**2.0+     &
      384.0*alpha**3.0*beta**2.0+384.0*alpha**2.0*beta**3.0+48.0*alpha**4.0*            &
      beta**2.0+144.0*alpha**2.0*beta**4.0-1056.0*alpha*beta**2.0+1008.0*d(0+1)         &
      *beta**2.0-96.0*alpha**2.0*beta+96.0*alpha**2.0*beta**2.0-48.0*c(0+1)*alpha**2.0  &
      -192.0*alpha*beta+192.0*c(0+1)*alpha**2.0*beta**2.0-192.0*d(0+1)*alpha**2.0*      &
      beta**2.0-768.0*d(0+1)*alpha*beta**2.0+12.0*c(0+1)-228.0*d(0+1)+64.0*beta**6.0    &
      +384.0*beta**5.0-96.0*d(0+1)**2.0-20.0*beta**4.0-1056.0*beta**3.0-96.0*           &
      alpha**3.0-12.0*alpha**4.0+240.0*alpha+240.0*beta-33.0*alpha**2.0+529.0*          &
      beta**2.0)/Dinf**2.0,                                                             &
      cmplx(-0.03125,0.0)*beta*d(0+1)-0.03125*alpha* d(0+1)-0.0625*beta**2*d(0+1)**2    &
      +0.125*d(0+1)*beta**3-0.125*alpha*beta**3-0.0078125*d(0+1)*alpha**2+(1.0/128)*    &
      c(0+1)*beta**2-0.0625*alpha**3*beta**2-0.0625*alpha**2*beta**3-0.0078125*alpha**4 &
      *beta**2-(1.0/128)*alpha**2*beta**4+0.015625*alpha*beta**2-0.0078125*d(0+1)*      &
      beta**2+0.015625*alpha**2* beta-(7.0/128)*alpha**2*beta**2+0.0078125*c(0+1)*      &
      alpha**2+0.03125*alpha*beta-0.03125*c(0+1)*alpha**2*beta**2+(1/32)*d(0+1)*        &
      alpha**2*beta**2+0.125*d(0+1)*alpha*beta**2-0.001953125*c(0+1)+0.001953125*       &
      d(0+1)+0.000244140625+(1/64)*d(0+1)**2-0.060546875*beta**4+0.015625*beta**3       &
      +0.015625*alpha**3+0.001953125*alpha**4-0.00390625*alpha-(1.0/256)*               &
      beta+0.01416015625*alpha**2+0.01416015625*beta**2 /), (/ 2,2 /)))  
      
      if (debug .eqv. .true.) then
       print *, 'Uleft(:,:,3,1):'
       call printRows(cmplxMatrix=Uleft(:,:,3,1))
      end if                 

      Uleft(:,:,3,2) = transpose(reshape((/ cmplx((-1.0/192.0),0.0)*beta**6.0+        &
                        35.0/768.0*beta**4.0-(259.0/3072.0) *beta**2.0+(75.0/4096.0), &
                        ((-1.0/192.0)*beta**6.0+(35.0/768.0)*beta**4.0-(259.0/3072.0) &
                        *beta**2.0+(75.0/4096.0))*cmplx(0.0,1.0)*Dinf**2.0,           &
                        cmplx(0.0,(-1.0/12288.0))*(64.0*beta**6.0-560.0*beta**4.0     &
                        +1036.0*beta**2.0-225.0)/cmplx(Dinf,0.0)**2.0,                &
                        cmplx((1.0/192.0),0.0)*beta**6.0-(35.0/768.0)*beta**4.0       &
                        +(259.0/3072.0)*beta**2.0-(75.0/4096.0) /), (/ 2,2 /)))
       if (debug .eqv. .true.) then
        print *, 'Uleft(:,:,3,2):'
        call printRows(cmplxMatrix=Uleft(:,:,3,2))
      end if   
  end subroutine compUExplicit


  subroutine compU(alpha, beta, Dinf, c ,d, maxOrder, Uleft, Uright)
    real(kind = PRECISION), intent(in) :: alpha
    real(kind = PRECISION), intent(in) :: beta
    real(kind = PRECISION), intent(in) :: Dinf
    complex(kind = PRECISION), intent(in), dimension(:) :: c
    complex(kind = PRECISION), intent(in), dimension(:) :: d
    integer, intent(in) :: maxOrder
    complex(kind = PRECISION), intent(out), dimension(:,:,:,:) :: Uleft
    complex(kind = PRECISION), intent(out), dimension(:,:,:,:) :: Uright
    complex(kind = PRECISION), dimension(:,:,:,:), allocatable :: Wr
    complex(kind = PRECISION), dimension(:,:,:,:), allocatable :: Wl
    integer :: mo, i,j,k,l,m,n

    mo = nint((maxOrder - 2.0) / 2)

    allocate(Wr(2, 2, maxOrder-1, mo+1 )); Wr = 0.0;
    allocate(Wl(2, 2, maxOrder-1, mo+1 )); Wl = 0.0;

    Wr = compWV(alpha,beta,Dinf,c,maxOrder,1,1)
    Wl = compWV(beta,alpha,Dinf,d,maxOrder,-1,1)

    do k = 1,(maxOrder-1)
      do m = 1,nint(k/2.0)
        Uright(:,:,k,m) = Wr(:,:,k,nint(k/2.0)+1-m);
        Uleft(:,:,k,m)  = Wl(:,:,k,nint(k/2.0)+1-m);
          do j=1,(k-1)
            do l = max(m-ceiling(j/2.0),1),ceiling((k-j)/2.0)
              Uright(:,:,k,m) = Uright(:,:,k,m) + &
              Uright(:,:,k-j,l)*Wr(:,:,j,nint(j/2.0)+1+l-m);
              Uleft(:,:,k,m) = Uleft(:,:,k,m) + &
              Uleft(:,:,k-j,l)*Wl(:,:,j,nint(j/2.0)+1+l-m);
            end do
          end do
          do j=1,(k-1)
            do n = 0,(ceiling(j/2.0)-m)
              do i = 1,ceiling((k-j)/2.0)
                Uright(:,:,k,m) = Uright(:,:,k,m) + &
                poch(1.0-i-n,n)/2**(i+n)/factorial(n)* &
                Uleft(:,:,k-j,i)*Wr(:,:,j,nint(j/2.0)+1-n-m)
                Uleft(:,:,k,m) = Uleft(:,:,k,m) + &
                poch(1.0-i-n,n)/(-2)**(i+n)/factorial(n)* &
                Uright(:,:,k-j,i)*Wl(:,:,j,nint(j/2.0)+1-n-m)
              end do
            end do
          end do
        end do
      end do
    deallocate(Wr); deallocate(Wl);
  end subroutine compU


  function compWV (q, t, Dinf, cd, maxOrder, r, isW) result(WV)
    real(kind = PRECISION), intent(in) :: q,t, Dinf
    complex(kind = PRECISION), intent(in), dimension(:) :: cd
    integer, intent(in) :: maxOrder, r , isW
    integer :: mo, n, i, j, k
    complex(kind = PRECISION) :: tmp
    complex(kind = PRECISION), dimension(:,:), allocatable :: g, st, stp, stm
    complex(kind = PRECISION), dimension(:), allocatable ::  &
                  f, OmOdd, OmEven, XiOdd, XiEven, ThOdd, ThEven
    complex(kind = PRECISION), dimension(:), allocatable ::  &
                        OmO, OmE, XiO, XiE, ThO, ThE
    complex(kind = PRECISION), dimension(:,:,:),   allocatable :: Ts
    complex(kind = PRECISION), dimension(:,:,:,:), allocatable :: WV
    complex(kind = PRECISION) :: a
    complex(kind = PRECISION) :: b
    logical, parameter :: debug = .false.

    if (debug .eqv. .true.) then
      print *, '---------compWV-----------------', r
      print *, 'cd: ', cd
    end if

    mo = nint((maxOrder - 2.0) / 2)

    allocate(f(mo+1))
    f = 0
    do n = 0,mo
      f(n+1) = poch(0.5, n) / (-r*2)**n / (1 + 2*n) / factorial(n)
    end do

    allocate(g(maxOrder-1,mo+1))
    g = cmplx(0.0,0.0)
    g(1,1) = cmplx(1.0,0.0)
    do n = 1,mo
    g(1,n+1) = dot_product(-g(1,1:(n+1)), flipVec(f(1:(n+1))))
    end do

    allocate(st(maxOrder-1,mo+1))
    allocate(stp(maxOrder-1,mo+1))
    allocate(stm(maxOrder-1,mo+1))
    do n=0,mo
      tmp = 0
      do j = 0,n
        tmp = tmp + binom(0.5, j) / (r*2.0)**j * cd(n-j+1)
      end do
      st(1,n+1)  = cmplx(r,0.0)*tmp + f(n+1)*(q+t)
      stp(1,n+1) = cmplx(r,0.0)*tmp + f(n+1)*(q+t+1)
      stm(1,n+1) = cmplx(r,0.0)*tmp + f(n+1)*(q+t-1)
    end do

    if (debug .eqv. .true.) then
      print *, 'st: ', st
      print *, 'stp: ', stp
      print *, 'stm: ', stm
    end if

    do i = 2,max(mo*2+1, (maxOrder-1))
      do n=0,mo
        g(i,n+1) = sum(g(i-1,1:(n+1))*flipVec(g(1,1:(n+1))) )
        st(i,n+1) = sum(st(i-1,1:(n+1))*flipVec(st(1,1:(n+1))) )
        stp(i,n+1) = sum(stp(i-1,1:(n+1))*flipVec(stp(1,1:(n+1))) )
        stm(i,n+1) = sum(stm(i-1,1:(n+1))*flipVec(stm(1,1:(n+1))) )
      end do
    end do

    allocate(OmOdd(mo+1));  OmOdd  = 0.0;
    allocate(OmEven(mo+1)); OmEven = 0.0;
    allocate(XiOdd(mo+1));  XiOdd  = 0.0;
    allocate(XiEven(mo+1)); XiEven = 0.0;
    allocate(ThOdd(mo+1));  ThOdd  = 0.0;
    allocate(ThEven(mo+1)); ThEven = 0.0;
    do n = 0,mo
      if ( (mod(maxOrder,2) .eq. 0) .or. (n .ne. mo) ) then
        OmEven(n+1) = OmEven(n+1) + st(1,n+1);
        XiEven(n+1) = XiEven(n+1) + stp(1,n+1);
        ThEven(n+1) = ThEven(n+1) + stm(1,n+1);
        do j = 1,n
            OmOdd(n+1) = OmOdd(n+1) + (r*2)**j*st(2*j,n-j+1)/factorial(2*j)
            XiOdd(n+1) = XiOdd(n+1) + (r*2)**j*stp(2*j,n-j+1)/factorial(2*j)
            ThOdd(n+1) = ThOdd(n+1) + (r*2)**j*stm(2*j,n-j+1)/factorial(2*j)
            OmEven(n+1) = OmEven(n+1) + (r*2)**j*st(2*j+1,n-j+1)/factorial(2*j+1)
            XiEven(n+1) = XiEven(n+1) + (r*2)**j*stp(2*j+1,n-j+1)/factorial(2*j+1)
            ThEven(n+1) = ThEven(n+1) + (r*2)**j*stm(2*j+1,n-j+1)/factorial(2*j+1)
        end do
      else !This splitting is needed to avoid wrong OmEven(n+1) etc when
           !n=mo and maxOrder is odd
        do j = 1,n
            OmOdd(n+1) = OmOdd(n+1) + (r*2)**j*st(2*j,n-j+1)/factorial(2*j)
            XiOdd(n+1) = XiOdd(n+1) + (r*2)**j*stp(2*j,n-j+1)/factorial(2*j)
            ThOdd(n+1) = ThOdd(n+1) + (r*2)**j*stm(2*j,n-j+1)/factorial(2*j)
        end do
      end if
    end do
    OmOdd(1) = cmplx(1.0,0.0)
    XiOdd(1) = cmplx(1.0,0.0) ! Overwrite because else wrong
    ThOdd(1) = cmplx(1.0,0.0)

    if (debug .eqv. .true.) then
      print *, 'OmEven', OmEven
      print *, 'XiEven', XiEven
      print *, 'ThEven', ThEven
    end if

    allocate(OmO(mo+1));  OmO = 0.0;
    allocate(OmE(mo+1));  OmE = 0.0;
    allocate(XiO(mo+1));  XiO = 0.0;
    allocate(XiE(mo+1));  XiE = 0.0;
    allocate(ThO(mo+1));  ThO = 0.0;
    allocate(ThE(mo+1));  ThE = 0.0;

    do n = 0,mo
      do j = 0,n
        OmO(n+1) = OmO(n+1) + binom(-0.5,j)/(r*2.0)**j*OmOdd(n-j+1)/sqrt(r*cmplx(2.0,0.0))
        XiO(n+1) = XiO(n+1) + binom(-0.5,j)/(r*2.0)**j*XiOdd(n-j+1)/sqrt(r*cmplx(2.0,0.0))
        ThO(n+1) = ThO(n+1) + binom(-0.5,j)/(r*2.0)**j*ThOdd(n-j+1)/sqrt(r*cmplx(2.0,0.0))
        if ( (mod(maxOrder,2) .eq. 0) .or. (n .ne. mo) ) then
            OmE(n+1) = OmE(n+1) + binom(-0.5,j)/(r*2.0)**j*OmEven(n-j+1)
            XiE(n+1) = XiE(n+1) + binom(-0.5,j)/(r*2.0)**j*XiEven(n-j+1)
            ThE(n+1) = ThE(n+1) + binom(-0.5,j)/(r*2.0)**j*ThEven(n-j+1)
        end if
      end do
    end do

    if (debug .eqv. .true.) then
      print *, 'OmO:', OmO
      print *, 'OmE', OmE
      print *, 'XiO', XiO
      print *, 'XiE', XiE
      print *, 'ThO', ThO
      print *, 'ThE', ThE
    end if

    allocate(Ts(2,2,mo+1)); Ts = 0.0;
    allocate(WV(2, 2, maxOrder-1, mo+1 )); WV = 0.0;

    do k = 1,(maxOrder-1)
      a = (q**2.0 + k/2.0 - 0.25)/k
      b = cmplx(0.0,-r*(k-0.5))
      mo = nint( (maxOrder-1-k)/2.0)-1 + nint(k/2.0)
      if (mod(k,2) .eq. 1) then
        do n = 0,mo
          Ts(1,1,n+1) = -a*binom(n-1.5,n)*(-r*2.0)**(-n)*(2.0*n+1) &
                        *sqrt(r/cmplx(2.0,0.0)) + cmplx(0.0,1.0)*b*OmO(n+1)
          Ts(1,2,n+1) = Dinf**2.0*(cmplx(0.0,1.0)*a*binom(-0.5,n)* &
                        (r*2.0)**(-n)/sqrt(r*cmplx(2.0,0.0)) +r*b*XiO(n+1))
          Ts(2,1,n+1) = (cmplx(0.0,1.0)*a*binom(-0.5,n)*(r*2.0)**(-n)/ &
                        sqrt(r*cmplx(2.0,0.0)) +r*b*ThO(n+1) )/Dinf**2.0
          Ts(2,2,n+1) = (a*binom(n-1.5,n)*(-r*2.0)**(-n)*(2.0*n+1.0)* &
                        sqrt(r/cmplx(2.0,0.0)) - cmplx(0.0,1.0)*b*OmO(n+1) )

        !print *, 'go frome here'
        !print *, 'n:', n
        !print *, 'brac:',    brac(k-1,q)/(r*2)**(3*k/2.0)
        !print *, 'g:',        g(k,1:(n+1))
        !print *, 'reshape:', reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /))
        !print *, 'repmat:',  repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),1)
        !print *, 'shape(repmat)', shape(repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),1))
        !print *, 'sum:', sum(repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),1),3)
        !print *, 'sum2:', sum(repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),1) * flipTensor(Ts),3)
        !print *, '-----------------------------------------------------'
        !print *, 'sum:',     sum(repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),2),3)
        WV(:,:,k,n+1) = brac(k-1,q)/(r*cmplx(2.0,0.0))**(3*k/cmplx(2.0,0.0)) &
                      * sum(repmat(reshape(g(k,1:(n+1) ),(/ 1,1,n+1 /)),1) &
                      * flipTensor(Ts),3)
        end do
     else
      !TODO: Translate later&
      !  do n = 0,mo
      !      Ts(:,:,n+1) = [(a*(n==0) +1i*r*b*OmE(n+1)), Dinf**2*b*XiE(n+1) ;  &
			!	b*ThE(n+1)/Dinf**2, (a*(n==0) -1i*b*r*OmE(n+1) )];
      !      WV(:,:,k,n+1) = brac(k-1)/(r*2)**(3*k/2)*sum(repmat(reshape(&
      !          g(k,1:(n+1) ),[1,1,n+1]),[2,2,1]).*Ts(:,:,(n+1):-1:1),3);
      !      if (isW .eq. 1) then ! The extra matrix of the theorem:
      !          WV(:,:,k,n+1) = WV(:,:,k,n+1) + brac(k-1)/(r*2)**(3*k/2)*( &
      !              -(4*q**2+2*k-1)*g(k,n+1)/2/k*eye(2) );
      !      end if
      !  end do
     end if
    end do

    if (debug .eqv. .true.) then
      print *, '----------Ts----------------'
      print *, Ts
      print *, '------WV------', shape(WV)
      print *, WV(:,:,1,1)
    end if

    deallocate(g);  deallocate(f);
    deallocate(st); deallocate(stp);  deallocate(stm);
    deallocate(OmOdd);  deallocate(OmEven); deallocate(XiOdd);
    deallocate(XiEven); deallocate(ThOdd);  deallocate(ThEven);
    deallocate(OmO);  deallocate(OmE); deallocate(XiO);
    deallocate(XiE); deallocate(ThO);  deallocate(ThE);
    deallocate(Ts);
  end function compWV

  !implement the rising pochhammer symbol. (briefly tested seems ok)
  function poch(x,n) result(res)
    real(kind = PRECISION), intent(in) :: x
    integer, intent(in) :: n
    real(kind = PRECISION) :: res
    integer :: i
    res = product(x + (/(i,i=0,(n-1))/))
  end function poch

  !implementation of the binomial coefficient.
  function binom(x,n) result(res)
    real(kind = PRECISION), intent(in) :: x
    integer, intent(in) :: n
    real(kind = PRECISION) :: res
    integer :: i
    res = product((/ ((x-i),i=0,n-1)  /)) / product((/ (i,i=1,n)  /))
  end function binom


  function brac(k,q) result(res)
    real(kind = PRECISION), intent(in) :: q
    integer, intent(in) :: k
    integer :: i
    real(kind = PRECISION) :: res

    if (k == 0) then
      res = 1
    else
      res = product(4*q**2 - (2*(/( i,i=1,k  )/) - 1)**2) /(2**(2*k)*(factorial(k) ))
    end if
  end function brac

  !An implementation for the factorial with overflow warning.
  function factorial(k) result(res)
    integer :: k
    real(kind = PRECISION) :: res
    integer :: i

    if (k > 170) then
      print *, ''//achar(27)//'[31m Error. Factorial overflow.', ''//achar(27)//'[0m'
    end if
    res = 1.0
    do i = 1,k
      res = res * i
    end do
  end function factorial

end module AsyModule
