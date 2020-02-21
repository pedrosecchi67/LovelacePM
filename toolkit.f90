!FORTRAN 90 subroutines to optimize mathematical operations, currently including:
! 1- AIC generation
! 2- Calculation of local self-induced velocity based on adjacent lines
! 3- Add wake influence to panel-wise AIC matrix

subroutine aicm_lines_gen(npan, nlin, lines, colpoints, aicm)
    integer, intent(IN) :: npan, nlin
    real(8), intent(IN) :: lines(1:nlin, 1:3, 1:2), colpoints(1:npan, 1:3)
    real(8), intent(OUT) :: aicm(1:3, 1:npan, 1:nlin)

    integer :: i, j
    real(8) :: a(3), b(3), na, nb

    do i=1, npan
        do j=1, nlin
            a=lines(j, 1:3, 1)-colpoints(i, 1:3)
            b=lines(j, 1:3, 2)-colpoints(i, 1:3)
            na=norm2(a)
            nb=norm2(b)
            aicm(1:3, i, j)=(((/a(2)*b(3)-a(3)*b(2), &
            a(3)*b(1)-a(1)*b(3), &
            a(1)*b(2)-a(2)*b(1)/)*(1.0/na+1.0/nb))/&
            (na*nb+dot_product(a, b)))/12.5663706
        end do
    end do
end subroutine aicm_lines_gen

subroutine self_influence(nlin, nloc, lines, solution, S, nvec, loclines, haswake, vdv)
    integer, intent(IN) :: nlin, nloc
    real(8), intent(IN) :: lines(1:nlin, 1:3, 1:2), solution(1:nlin), S, nvec(1:3)
    integer, intent(IN) :: loclines(1:nloc)
    logical, intent(IN) :: haswake
    real(8), intent(OUT) :: vdv(3)

    integer :: i
    real(8) :: Gamma(3)

    vdv=(/0.0, 0.0, 0.0/)

    do i=1, nloc
        Gamma=solution(loclines(i))*(lines(loclines(i), 1:3, 2)-lines(loclines(i), 1:3, 1))
        vdv=vdv+(/Gamma(2)*nvec(3)-Gamma(3)*nvec(2), Gamma(3)*nvec(1)-Gamma(1)*nvec(3), Gamma(1)*nvec(2)-Gamma(2)*nvec(1)/)
    end do
    
    if(haswake) then
        vdv=vdv/((nloc+1)*S)
    else
        vdv=vdv/(nloc*S)
    end if
end subroutine self_influence