!FORTRAN 90 subroutines to optimize mathematical operations, currently including:
! 1- AIC generation

subroutine gen_aicm(npan, nwake, coords, colpoints, addto, aicm)!x, aicmy, aicmz)
    integer, intent(IN) :: npan, nwake
    real(8), intent(IN) :: coords(1:npan+nwake, 1:3, 1:4), colpoints(1:npan, 1:3)
    integer, intent(IN) :: addto(1:nwake, 1:2)
    real(8), intent(OUT) :: aicm(1:3, 1:npan, 1:npan)!x(1:npan, 1:npan), aicmy(1:npan, 1:npan), aicmz(1:npan, 1:npan)

    integer :: i, j
    real(8) :: p1(3), p2(3), np1, np2, vbar(3)

    do i=1, npan
        do j=1, npan+nwake
            vbar=(/0.0, 0.0, 0.0/)

            !segment 1
            p1=coords(j, 1:3, 1)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 2)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 2
            p1=coords(j, 1:3, 2)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 3)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 3
            p1=coords(j, 1:3, 3)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 4)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 4
            p1=coords(j, 1:3, 4)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 1)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            if(j<=npan) then
                aicm(1:3, i, j)=vbar
            else
                if(addto(j-npan, 1)/=0) then
                    aicm(1:3, i, addto(j-npan, 1))=aicm(1:3, i, addto(j-npan, 1))+vbar
                end if
                if(addto(j-npan, 2)/=0) then
                    aicm(1:3, i, addto(j-npan, 2))=aicm(1:3, i, addto(j-npan, 2))-vbar
                end if
            end if
        end do
    end do
end subroutine gen_aicm

subroutine gen_aicm_nowake(npan, coords, colpoints, aicm)!x, aicmy, aicmz)
    integer, intent(IN) :: npan
    real(8), intent(IN) :: coords(1:npan, 1:3, 1:4), colpoints(1:npan, 1:3)
    real(8), intent(OUT) :: aicm(1:3, 1:npan, 1:npan)!x(1:npan, 1:npan), aicmy(1:npan, 1:npan), aicmz(1:npan, 1:npan)

    integer :: i, j
    real(8) :: p1(3), p2(3), np1, np2, vbar(3)

    do i=1, npan
        do j=1, npan
            vbar=(/0.0, 0.0, 0.0/)

            !segment 1
            p1=coords(j, 1:3, 1)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 2)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 2
            p1=coords(j, 1:3, 2)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 3)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 3
            p1=coords(j, 1:3, 3)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 4)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            !segment 4
            p1=coords(j, 1:3, 4)-colpoints(i, 1:3)
            p2=coords(j, 1:3, 1)-colpoints(i, 1:3)
            !p1
            if(norm2(p2-p1)>0.00001) then
                np1=norm2(p1)
                np2=norm2(p2)
                vbar=vbar+((1.0/np1+1.0/np2)*(/p1(2)*p2(3)-p1(3)*p2(2), &
                p1(3)*p2(1)-p1(1)*p2(3), p1(1)*p2(2)-p1(2)*p2(1)/)/(np1*np2+ &
                dot_product(p1, p2)))/12.5663706
            end if

            aicm(1:3, i, j)=vbar
        end do
    end do
end subroutine gen_aicm_nowake