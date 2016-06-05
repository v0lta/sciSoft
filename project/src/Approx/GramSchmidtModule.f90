module GramSchmidtModule
  use UtilModule
  use FunctionModule
  use IntegrationModule
  implicit none
  public

contains

  !compute c using least squares.
  function compC(A, f, row, col) result(c)
    integer :: row, col
    type(functionType) :: f
    real(kind = PRECISION), dimension(:,:)  :: A
    real(kind = PRECISION), dimension(row)  :: b
    real(kind = PRECISION), dimension(col) :: c
    integer :: rank, info, ldb
    integer :: lwork
    real :: rcond = 0.001
    real(kind = PRECISION), dimension(:), allocatable :: S
    real(kind = PRECISION), dimension(:), allocatable :: work

    ldb = max(row,col)
    lwork = 4 * row + 10*max(2*row,col)
    allocate(work(lwork))
    allocate(S(min(row,col)))
    b = f%functionValues
    !c = pinv(P) * f;
    !SUBROUTINE SGELSS( M,  N, NRHS,A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )
    call        DGELSS( row,col, 1,   A, row,   b, ldb, S, rcond, rank, work, lwork, info)
    !print *, 'lapack says:', info
    c = b(1:col)
    deallocate(work)
    deallocate(S)
  end function compC

  !Orthogonalize using Gram-Schmidt.
  subroutine gramSchmidt(N,gAlpha,gBeta,wAlpha,wBeta,polyMat,meth)
        integer, intent(in) :: N, meth
        integer :: columns, i, j, k
        real(kind = PRECISION), intent(in) :: gAlpha, gBeta, wAlpha, wBeta
        real(kind = PRECISION), intent(inout), dimension(:,:) :: polyMat
        real(kind = PRECISION), dimension(N) :: x, orthFunVals, normFunVals
        real(kind = PRECISION) ::  eps, intRes, parPart
        procedure (funPointer), pointer :: fpOrth => null()
        procedure (funPointer), pointer :: fpNorm => null()

        !nagfor, ifort etc really dont like the f2008 standard yet. Thus initialization
        !cannot happen at declaration.
        fpOrth => orthFun
        fpNorm => normFun

        !Public variables allow for interconnection with the quadpack inter-
        !face.
        pAlpha = wAlpha
        pBeta = wBeta
        allocate(pol1(size(polyMat,1)))
        allocate(pol2(size(polyMat,1)))

        columns = size(polyMat,2)
        call linspace(-1.0,1.0,x)

        !Transform to remove the extrema and be able to handle negative
        !alpha and beta.
        eps = sqrt(epsilon(x))
        x = (x - minval(x))/(maxval(x) - minval(x)) * (2 - 2*eps) + eps - 1

        !Orthogonalize
        do i = 1,columns
          do j=1,i
            if (i .ne. j) then
              !Integrate the orthogonalty function.
              !Set the public variables.
              pol1 = polyMat(:,i)
              pol2 = polyMat(:,j)
              if (meth .eq. 1) then
                !use gaussian quadrature,
                do k = 1,N
                  orthFunVals(k) = orthFun(x(k))
                  normFunVals(k) = normFun(x(k))
                end do
                parPart = integral(gAlpha,gBeta,N,orthFunVals) / &
                          integral(gAlpha,gBeta,N,normFunVals)
              elseif (meth .eq. 2) then
                !use quadpack.
                parPart = quadIntegral(fpOrth, minval(x), maxval(x)) / &
                          quadIntegral(fpNorm, minval(x), maxval(x))
              else
                !use the trapeziodal rule.
                parPart = trapzPoint(fpOrth, minval(x), maxval(x), N) / &
                          trapzPoint(fpNorm, minval(x), maxval(x), N)
              end if
              polyMat(:,i) = polyMat(:,i) - parPart * polyMat(:,j)
            end if
          end do
          !normalize:
          !update the public row variable
          pol1 = polyMat(:,i)
          !use gaussian quadrature-integration.
          if (meth .eq. 1) then
            !Evaluate the norm-Functoin.
            do k = 1,N
                normFunVals(k) = normFun(x(k))
            end do
            intRes = integral(gAlpha,gBeta,N,normFunVals)
          elseif (meth .eq. 2) then
            !use quadpack.
            intRes = quadIntegral(fpNorm, minval(x), maxval(x))
          else
            !use the trapeziodal rule.
            intRes = trapzPoint(fpNorm, minval(x), maxval(x), N)
          end if
            polyMat(:,i) = polyMat(:,i) / sqrt(intRes)
        end do

        deallocate(pol1)
        deallocate(pol2)
  end subroutine gramSchmidt

  !The normalization and orthogonalty functions required by GM.
  function orthFun(x) result(res)
      real(kind = PRECISION), intent(in) :: x
      real(kind = PRECISION) :: res

      res = omega(x,pAlpha,pBeta) * polyval(pol1,x) * polyval(pol2,x)
  end function orthFun

  function normFun(x) result(res)
    real(kind = PRECISION), intent(in) :: x
    real(kind = PRECISION) :: res

    res = omega(x,pAlpha,pBeta) * polyval(pol1,x) ** 2
  end function

  !The weight function.
  function omega(x,alpha,beta) result(res)
    real(kind = PRECISION) :: x, alpha, beta, res
    !!for -1 this thing blows up!!!
    if (x == -1) then
      res = 0
    else
      res = (1 - x)**alpha * (1 + x)**beta
    end if
  end function omega

  !Transform the coordinates to -1,1
  subroutine transformX(x,des, transA,transB)
    real(kind = PRECISION), dimension(:), intent(inout) :: x, des
    real(kind = PRECISION), intent(out) :: transA,transB

    transB = 2/(maxval(x) - minval(x))
    transA = (-1/transB) - minval(x)
    x = (x + transA)*transB
    des = (des + transA)*transB
  end subroutine transformX

  !Undo the transformation to -1,1
  subroutine reverseTransformX(x,des, transA,transB)
    real(kind = PRECISION), dimension(:), intent(inout) :: x, des
    real(kind = PRECISION), intent(in) :: transA,transB

    x = x/transB - transA
    des = des/transB - transA
  end subroutine reverseTransformX

  !Compute the matrix, which contains the evaluated polynomials in the
  !desired range for x.
  subroutine evalMatrix(polyMat,x,M,evalMat)
    real(kind = PRECISION), dimension(:,:), intent(inout) :: polyMat
    real(kind = PRECISION), dimension(:,:), intent(inout) :: evalMat
    real(kind = PRECISION), dimension(:), intent(in) :: x
    integer, intent(in) :: M
    integer :: i, k, length

    length = size(x)
    do i = 1,M
      do k= 1,length
        evalMat(k,i) = polyval(polyMat(:,i), x(k))
      end do
    end do
  end subroutine evalMatrix
end module
